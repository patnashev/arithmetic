# State serialization — deep dive

Every checkpoint, recovery point, proof point, certificate, and progress file a consumer application writes goes through one small layer: `Writer`/`Reader`/`File` in `file.{h,cpp}`. A task's `TaskState` subclass describes *what* to serialize (its fields); `File` wraps that in a fixed binary header, hashes it, and writes it atomically. This is the **on-disk contract**: any change to the header layout, the `TYPE` discriminators, or a `TaskState`'s field order breaks resumability for in-flight work and for files written by older builds.

Three facts to anchor on before reading the code:

1. **Files are typed and fingerprinted.** Every file starts with a magic number, an `appid`, a one-byte `TYPE`, a version byte, and (usually) a 32-bit fingerprint. `File::read(TaskState&)` refuses to load a buffer whose `TYPE` doesn't match the state being read into — so a `.ckpt` from one test class can't be misread as another's.
2. **`TYPE` is a per-application, hand-assigned namespace.** The registry of values is owned by the consumer: each application keeps its own set of `static const` TYPE constants on its `TaskState` subclasses. PRST's registry is documented in `checkpoints.md` (in the patnashev/prst repo); Prefactor uses its own checkpoints — they're completely different. Reusing or renumbering a value within an application silently aliases two on-disk formats.
3. **Writes are crash-safe; reads are hash-checked.** `File::commit_writer` writes to `<name>.new` (on Windows via `FILE_FLAG_WRITE_THROUGH` so the data is flushed to disk; on POSIX via plain `fopen`+`fwrite`+`fclose` with **no** `fsync`, so durability isn't guaranteed before the rename), then removes the old file and renames `.new` over the original, and drops a `<name>.md5` sidecar. `read_buffer` verifies that sidecar and discards a corrupt buffer (treated as "no checkpoint" → restart from scratch).

Source files:
- `file.h`, `file.cpp` (`Writer`, `Reader`, `TextReader`, `File` + `FileEmpty`/`FilePacked`)
- `md5.{h,c}` (the MD5 used for the sidecar and `unique_fingerprint`)
- `task.{h,cpp}` (the `TaskState` base class — §4)

Prereqs / companions: `task-lifecycle.md` (the `Task` runtime; the `on_state` cadence that decides *when* a checkpoint is written), `checkpoints.md` (in the patnashev/prst repo) (PRST's TYPE registry, its `.ckpt`/`.rcpt` strong-check files, and the LLR2 compatibility layer), `proof-system.md` (in the patnashev/prst repo) (the proof records and the `.pack` `FilePacked` container), `boinc-and-net.md` (in the patnashev/prst repo) (`NetFile` ships these exact buffers over HTTP).

## 1. The on-disk layout

Every file `File::write` produces has this shape:

```
offset  size  field
  0      4    MAGIC_NUM = 0x9f2b3cd4          (uint32, little-endian)
  4      1    appid                            (PRST = 4; default 1; LLR2 = 2)
  5      1    format_version                   (always 0 today — the "(0 << 8)")
  6      1    TYPE                             (the TaskState discriminator)
  7      1    version                          (TaskState::version(), always 0 today)
  8      4    fingerprint                      (uint32 — present only if fingerprint != 0)
 8/12   ...   body                             (the TaskState record)
```

Bytes 4–7 are a single `uint32` packed as `appid | (format_version << 8) | (type << 16) | (version << 24)` (`file.cpp:205`). The fingerprint word is conditional: a file constructed with fingerprint `0` (e.g. progress `.param` files) omits it and the body starts at offset 8; otherwise the body starts at offset 12.

The **body** is whatever the `TaskState` subclass writes, and it always begins with the base class's 4-byte `iteration` (`task.cpp:21-24`), followed by the subclass's own fields in declaration order.

A separate **`<filename>.md5`** sidecar holds the 32-char hex MD5 of the file body (when `File::hash` is true, the default).

## 2. The write/read round trip

**Writing** (`file.cpp:355-369`, header at `:201-209`):

```cpp
Writer* File::get_writer(char type, char version)
{
    Writer* writer = get_writer();                  // reuses the file's buffer, cleared
    writer->write(MAGIC_NUM);
    writer->write(appid + (0 << 8) + ((unsigned char)type << 16) + ((unsigned char)version << 24));
    if (_fingerprint != 0)
        writer->write(_fingerprint);
    return writer;
}

bool File::write(TaskState& state)
{
    try {
        std::unique_ptr<Writer> writer(get_writer(state.type(), state.version()));
        state.write(*writer);                        // body: iteration + subclass fields
        commit_writer(*writer);                      // atomic .new → rename, + .md5 sidecar
    } catch (const std::exception&) { return false; }
    state.set_written();
    return true;
}
```

**Reading** (`file.cpp:259-275, 342-353`):

```cpp
Reader* File::get_reader()
{
    read_buffer();                                   // load file + verify .md5
    if (_buffer.size() < 8) return nullptr;
    if (*(uint32_t*)_buffer.data() != MAGIC_NUM) return nullptr;
    if (_buffer[4] != appid) return nullptr;                       // appid must match
    if (_fingerprint != 0 && _buffer.size() < 12) return nullptr;
    if (_fingerprint != 0 && *(uint32_t*)(_buffer.data() + 8) != _fingerprint) return nullptr;  // fingerprint must match
    return new Reader(_buffer[5], _buffer[6], _buffer[7], _buffer.data(), (int)_buffer.size(),
                      _fingerprint != 0 ? 12 : 8);   // pos = start of body
}

bool File::read(TaskState& state)
{
    std::unique_ptr<Reader> reader(get_reader());
    if (!reader) return false;
    if (reader->type() != state.type()) return false;   // TYPE must match the target state
    if (!state.read(*reader)) return false;              // truncation/format error → false
    state.set_written();
    return true;
}
```

Every gate that fails returns `nullptr`/`false`, which the caller treats as "no usable checkpoint" → the test starts (or restarts) from the beginning. There is no partial recovery: a file that fails magic, appid, fingerprint, TYPE, or a truncated body read is simply ignored. That's the safe default — a stale or foreign file never corrupts a run.

## 3. Field & method reference

**`Writer`** (`file.h:13-40`, `file.cpp:17-80`) — append-only byte buffer:

| Method | Encoding |
|---|---|
| `write(T)` (arithmetic) | raw `sizeof(T)` bytes, little-endian (template, `file.h:23-24`) |
| `write(const std::string&)` | `uint32` length + raw bytes (`file.cpp:22-26`) |
| `write(const Giant&)` | `int32` length **(negated if value < 0 — the sign carries the Giant's sign)** + `len*4` limb bytes (`file.cpp:28-35`) |
| `write(const SerializedGWNum&)` | `uint32` word count + `count*4` bytes (`file.cpp:37-41`) |
| `write_text` / `write_textline` | raw text; `_textline` appends `\r\n` |
| `hash()` / `hash_str()` | MD5 of the buffer (16 bytes / 32 hex chars) |

**`Reader`** (`file.h:42-65`, `file.cpp:82-157`) — positional reader over the loaded buffer, carrying the parsed `type`/`version`. Each `read` bounds-checks against `_size` and returns `false` on underrun. The `Giant` reader reconstructs sign from the negated length (`:140-141`); the `SerializedGWNum` reader takes word count × 4 bytes.

**`File`** (`file.h:82-126`) — the unit of persistence:

| Member | Role |
|---|---|
| `read(TaskState&)` / `write(TaskState&)` | the typed round trip above. |
| `get_reader` / `get_writer(type,version)` | low-level header framing. |
| `get_textreader` / `write_text` | unframed text (used by `-batch` files, the `stderr` net log). |
| `read_buffer` / `commit_writer` / `free_buffer` / `clear` | the I/O verbs (overridden by subclasses). |
| `add_child(name, fp)` | a dotted-name sub-file (e.g. per-base `.ckpt.<a>`), inheriting `hash`. |
| `unique_fingerprint(fp, id)` | `static`: `MD5(fp-bytes ‖ id)`, first 4 bytes → `uint32` (`file.cpp:182-192`). |
| `hash` (bool) | whether to write/verify the `.md5` sidecar (default true). |
| `appid` (int) | the byte-4 application id; PRST sets `File::FILE_APPID = 4` (`prst/src/prst.cpp:53`). |

**The `File` hierarchy:**
- `File` — local disk (`read_buffer` = `fopen`+read+verify-md5; `commit_writer` = atomic write).
- `FileEmpty` (`file.h:128-136`) — a `/dev/null`: reads nothing, writes nothing. Used where a `File&` is required but no persistence is wanted.
- `FilePacked` (`file.h:143-158`, `file.cpp:385-509`) — backed by a `container::FileContainer` (the proof `.pack`); `get_writer` returns a streaming `FilePackedWriter`. See `proof-system.md`.
- Consumers can subclass `File` for foreign-format interop — PRST's `LLR2File` is the worked example (`checkpoints.md` in the patnashev/prst repo).

## 4. The `TaskState` record & the checkpoint lifecycle

`File::read` keys off `TYPE`, so the constant is the format's identity. The base class (`task.h:23-43`):

```cpp
class TaskState {
    char _type;           // TYPE constant, sets the on-disk record identity
    bool _written;        // whether this state has been flushed to _file already
    int  _iteration;      // current iteration count (the only field the base class persists)

    void set(int iteration);      // resets _written to false and stores _iteration
    virtual bool read(Reader&);   // base reads _iteration only
    virtual void write(Writer&);  // base writes _iteration only
    void set_written();           // mark state as flushed; Task::write_state skips re-writing
    char type();                  // returns _type (the on-disk discriminator)
    char version() { return 0; }  // per-TYPE version byte; reserved for future schema bumps
    bool is_written();             // returns _written
    int  iteration();              // returns _iteration
};

template<class State>
State* read_state(File* file);   // null-safe deserializer
```

Subclasses override `read()`/`write()` and add their own fields; the record is always the base class's `iteration (uint32)` followed by the subclass's fields in declaration order. Conventions:

- `read_state<T>(file)` is the canonical way to load a checkpoint: it returns a fresh `T` if the file decodes, or `nullptr` — which callers treat as "no checkpoint".
- The `_written` flag prevents `Task::write_state()` (`task.cpp:192`) from re-writing a state that's already on disk; `TaskState::set(...)` clears it.
- Subclasses' `set(...)` methods take the same `int iteration` first argument plus the state-specific payload (a `Giant`, a `SerializedGWNum`, …).
- `bool`s are written as `int` `1`/`0`, not a single byte.
- `version()` is `0` for every state today (`task.h:35`), so header byte 7 is always 0. The version byte exists precisely so a future field change can be made readable old-and-new without a TYPE bump — bump `version()` and branch in `read`.
- **A state need not have a subclass at all.** A bare `TaskState(TYPE)` persists only the iteration — a placeholder record whose presence and iteration number carry all the meaning. PRST's TYPE 5 works this way (see `checkpoints.md` in the patnashev/prst repo).

A checkpoint's lifecycle (see `task-lifecycle.md` for the cadence): `Task` calls `File::write(state)` on the `DISK_WRITE_TIME` interval and when `Task::write_state()` is called. `write_state()` is public (`task.h:77`) and can be called after `run()` to persist the result — there is no automatic write at completion; checkpoints are usually erased on completion (`File::clear`). On restart `File::read(state)` repopulates `state` and the test resumes from `state.iteration()`.

## 5. Integrity, atomicity, and fingerprints

**Atomic write** (`file.cpp:283-323`): `commit_writer` writes the buffer to `<name>.new` via `writeThrough` (Windows `FILE_FLAG_WRITE_THROUGH`, else `fopen`+`fwrite`), `remove`s the old file, then `rename`s `.new` over it. A crash mid-write leaves an intact original plus a stray `.new`; it never leaves a half-written checkpoint. The `.md5` sidecar is written after.

**MD5 verification** (`file.cpp:234-256`): `read_buffer` loads the file, then — if `hash` and a non-empty `<name>.md5` exists — recomputes the MD5 and **clears the buffer on mismatch**. A corrupt checkpoint is thus indistinguishable from a missing one: the test restarts rather than trusting bad state. `md5.{h,c}` provides `MD5Init/Update/Final` and the `md5_raw_input(out[33], …)` hex helper.

**`unique_fingerprint`** (`file.cpp:182-192`): `MD5(fingerprint-bytes ‖ unique_id)[0:4]` as a `uint32`. This is how per-base / per-point child files get distinct fingerprints (e.g. `Pocklington` opens `.ckpt.<a>` children, the proof system salts by point position) so two sub-runs can't read each other's checkpoints even under the same parent filename.

**Foreign-format interop** is done by subclassing `File` and munging the header bytes in `read_buffer`/`commit_writer` overrides. PRST's `LLR2File`/`LLR2NetFile` (bidirectional compatibility with LLR2's on-disk format) is the worked example — documented in `checkpoints.md` (in the patnashev/prst repo). Byte 4 is the **appid**, not a format version; that's what such a subclass rewrites.

## 6. Pitfalls

- **`TYPE` is a permanent on-disk contract.** Renumbering or reusing a value silently makes two record formats collide; `File::read`'s type-match (`file.cpp:347`) would accept the wrong body and mis-deserialize. Consult your application's registry before assigning (PRST's is in `checkpoints.md` in the patnashev/prst repo) and only ever *append* new TYPEs. Bumping a TYPE without a backward-compat read path orphans every in-flight checkpoint.
- **Field order is the format.** `Writer`/`Reader` are positional with no field tags. Inserting a field in the middle of a `TaskState::write`, or reordering, breaks every existing checkpoint for that TYPE. Append at the end and gate on `version()` if you must extend.
- **`Reader::read(double)` under-checks bounds.** It verifies only 4 bytes remain but reads and advances 8 (`file.cpp:109-116`) — a latent over-read if a double sits in the last 4–7 bytes. No current `TaskState` serializes a `double`, so it's dormant. A fix is planned; verify against the current `file.cpp` before relying on this note.
- **appid is part of the identity.** A file written with the default `appid 1` or by LLR2 (`2`) won't load under PRST (`appid 4`) without the `LLR2File` path — `get_reader` rejects on `_buffer[4] != appid` (`file.cpp:267`). If checkpoints "aren't resuming" across tools, check the appid byte first.
- **`hash = false` files skip integrity *and* the sidecar.** Progress `.param` files set `hash = false` (they're frequently rewritten); they get no `.md5` and no corruption check. Don't assume every PRST file is hash-protected.
- **A fingerprint-0 file has no fingerprint word.** The body offset shifts from 12 to 8. Code that pokes a fixed offset (like the offset-12 iteration in PRST's LLR2 munging — `checkpoints.md` in the patnashev/prst repo) implicitly assumes a fingerprinted file.

## 7. Quick reference

Header bytes: `[0:4]` magic `0x9f2b3cd4` · `[4]` appid (PRST 4) · `[5]` format_version (0) · `[6]` TYPE · `[7]` version (0) · `[8:12]` fingerprint (if ≠ 0) · then body.

| You want to… | Where / how |
|---|---|
| Add a new checkpointable state | new `TaskState` subclass, **append** an unused `TYPE` from your application's registry (PRST: `checkpoints.md`), implement `read`/`write` calling `TaskState::` first |
| Extend an existing record | append fields at the end + bump `version()` + branch in `read` (don't reorder) |
| Make a file skip hashing | set `File::hash = false` before writing |
| Give a child file a distinct id | `File::unique_fingerprint(parent_fp, "<salt>")` |
| Persist nothing (but satisfy a `File&`) | `FileEmpty` |
| Pack many files into one blob | `FilePacked` + `container::FileContainer` (proof `.pack`) |
| Interop with a foreign checkpoint format | subclass `File`, munge the header in `read_buffer`/`commit_writer` (PRST's LLR2 layer: `checkpoints.md`) |
| Persist a bare "I was here" marker | `TaskState(TYPE)` with no subclass — iteration-only record (§4) |
| Encode a signed bignum | `Writer::write(Giant)` — sign rides the length field |

## 8. Open questions / non-coverage

- **`SerializedGWNum` internals.** This doc treats it as "uint32 word count + words" on the wire (`file.cpp:37-41, 146-157`); what those words *mean* (the GWnum `gwserialize` representation) is a GWnum concern, covered at the binding level in `arithmetic-foundation.md` §5. Whether a serialized value is portable across FFT *lengths* for the same modulus (and thus whether a `SerializedGWNum`-only checkpoint survives a reliable-mode FFT bump) is an open question flagged there — do **not** assume the earlier wording that it's "only meaningful to a compatible FFT setup"; that was unverified.
- **`container::FileContainer` / `Packer`.** The `.pack` streaming format (`container.{h,cpp}`, ~1800 lines of JSON-framed binary with its own MD5/codec/recovery) is summarized in `proof-system.md` §6 as a contract; the implementation is still unowned by any doc.
- **The MD5 implementation** (`md5.c`) is stock RFC 1321; not reviewed here beyond its interface.
- **Endianness.** All multi-byte writes are raw host-order (`write((char*)&value, sizeof)`), so checkpoints are implicitly little-endian and not portable to a big-endian host. PRST targets x86/ARM-LE, so this is academic — but it's an unstated assumption, not a guarantee.
- **The `format_version` byte (offset 5)** is always 0 and never branched on; it's a reserved second escape hatch alongside the per-TYPE `version` byte. Its intended use vs. `version()` isn't documented in code.
