# Arithmetic Framework — Lay of the Land

A high-level orientation to this library for contributors new to the codebase. This repo
(`patnashev/arithmetic`) is the **shared arithmetic framework**: the bignum, arithmetic, checkpointing,
logging, and option-parsing primitives that higher-level tools reuse. Its primary consumer is **PRST**
(`patnashev/prst`, a primality-testing utility), which embeds this repo as a `framework/` submodule;
other Atnashev-family tools consume it the same way. For deeper dives into each subsystem see the
per-subsystem docs in this folder; the application that drives them is documented in the prst repo's
`docs/` (start from its `lay-of-the-land.md`).

## What this library is

A C++ framework for **working with very large numbers** and the machinery around long-running
computations on them. It splits into roughly three layers:

- **Arithmetic substrate** (`arithmetic/`) — `Giant` (a heap bignum, ABI-compatible with gwnum's C
  `giant`) with `Giants`/`GMP`/`GW` backends; `GWArithmetic`/`GWNum`/`GWState` (FFT-based modular
  arithmetic over George Woltman's GWnum); number-theory helpers in `integer.*` (gcd, inv, is_prime,
  phi, jacobi, kronecker, a prime sieve); and specialized arithmetics — `lucas.*`, `edwards.*`,
  `montgomery.*`, `poly.*`, `group.*`.
- **Task & I/O layer** (repo root) — the `Task`/`InputTask` base classes (restart, checkpoint cadence,
  abort handling) that every long computation subclasses; `Logging`/`SubLogging`/`Progress` for results,
  progress, and time accounting; `File` for fingerprinted checkpoint storage; `FileContainer` for the
  `.pack` bundle format; `md5` hashing.
- **Front-end helper** — the `Config` option-parser DSL that consumers use to turn `argv`/`.ini` into
  their own options structs.

The library owns no `main()`; it is always driven by a consumer.

## Repo layout

```
arithmetic/                          ← repository root (patnashev/arithmetic)
├── README.md
├── docs/                            ← this folder (framework docs + mult_*.pdf references)
├── arithmetic/                      ← the arithmetic substrate
│   ├── giant.cpp / giant.h          ← Giant bignum + Giants/GMP/GW backends
│   ├── arithmetic.cpp / .h          ← GWState, GWArithmetic, GWNum, SerializedGWNum
│   ├── integer.cpp / .h             ← number theory (gcd/inv/is_prime/phi/jacobi/…), prime sieve
│   ├── lucas.cpp / .h               ← Lucas-sequence arithmetic (V/U recurrences)
│   ├── edwards.cpp / .h             ← Edwards-curve arithmetic (ECM; dormant in PRST)
│   ├── montgomery.cpp / .h          ← Montgomery-curve arithmetic (ECM; dormant in PRST)
│   ├── poly.cpp / .h                ← polynomial multiplication (dormant in PRST)
│   ├── group.cpp / .h               ← curve/poly group-formula implementations (large)
│   ├── field.h / exception.h        ← field abstraction; exception types
│   └── test.cpp                     ← standalone arithmetic smoke test (own main())
├── inputnum.cpp / .h                ← number parser + form classification (KBNC/FACTORIAL/…)
├── logging.cpp / .h                 ← Logging / SubLogging / Progress
├── task.cpp / .h                    ← Task / InputTask base classes
├── file.cpp / .h                    ← File checkpoint abstraction (+ FilePacked bridge)
├── config.cpp / .h                  ← option-parser DSL
├── container.cpp / .h               ← FileContainer (.pack proof bundles)
├── md5.c / md5.h                    ← hashing
├── gwnum/                           ← prebuilt GWnum (Woltman's FFT-based bignum; third-party)
├── gmp/                             ← prebuilt GMP for Windows (third-party)
└── bow/                             ← Big Object Writer (BOINC-only state-serialization glue)
```

Mental model: **this library is "the shared bignum + arithmetic + checkpointing primitives every
computation reuses"; the consumer (e.g. PRST) decides "what computation do we run and how do we
orchestrate it?"** Because the surface here is shared across tools, it is a higher-cost change surface
than a consumer's own code — read the relevant deep-dive before extending an interface.

## The pieces, and where to read more

- **`Giant` / `GWState` / `GWArithmetic`** — the bignum heap model and the FFT math runtime
  (instruction set, FFT config, modulus, fingerprint, `known_factors` reduction), plus `SerializedGWNum`
  (the checkpoint bridge) and `ReliableGWArithmetic` round-off escalation. See `arithmetic-foundation.md`.
- **Curves & polynomials** — the full (`GroupArithmetic`, NAF) vs. differential
  (`DifferentialGroupArithmetic`, DAC ladder) split; the five concrete arithmetics; chain construction.
  Only `LucasV`/`LucasUV` is a live PRST consumer; the Edwards/Montgomery ECM and `PolyMult` paths are
  dormant in PRST. See `curves-and-polynomials.md`.
- **`Task` / `InputTask`** — the restart loop in `run()`, the checkpoint/progress cadence in
  `on_state()`, the error-correction handshake, and the `TaskState` base. Every long computation
  subclasses these. See `task-lifecycle.md`.
- **`Logging` / `SubLogging` / `Progress`** — results/factor/log files, the progress checkpoint, the
  parent walk in `update()`, and the staged time accounting (`_time_total`/`_time_stage`). A subtle hot
  zone. See `logging-and-progress.md`.
- **`InputNum`** — parse + form classification (`KBNC`/`FACTORIAL`/`PRIMORIAL`/`GENERIC`), the
  cofactor-divisor form, algebraic factoring, `mod`/fingerprinting, `is_half_factored`. See
  `inputnum-parsing.md`.
- **`Writer` / `Reader` / `File` serialization** — the on-disk format (magic + appid + `TYPE` + version
  + fingerprint header), atomic write + `.md5` sidecar, `unique_fingerprint`, and the appid scheme that
  lets consumers share the format. See `state-serialization.md`.
- **`Config` DSL** — the runtime `ConfigObject` tree vs. the CRTP `*Setup` builder; `delim` matching;
  the builder verbs; the `parse_args` walk and `.ini` parsing. See `config-dsl.md`.
- **`FileContainer` / `.pack`** — the bundle format (JSON-framed chunks, stream index, codecs, MD5
  recovery) and the `FilePacked` File-layer bridge. See `container-format.md`.

## Build

This library is consumed as a submodule and built into the consumer's binary; it has no standalone build
target beyond `arithmetic/test.cpp` (compiled and run by hand for smoke-testing). GWnum and GMP ship as
prebuilt static libraries under `gwnum/{linux64,mac64,win64}` and `gmp/win64`; a consumer that needs a
newer GWnum bumps the submodule pointer. Per-platform build wiring (solutions/makefiles) lives with the
consumer, not here.

## Companion deep-dives in this folder (framework subsystems)

- **`task-lifecycle.md`** — `TaskState` / `Task` / `InputTask`, the restart loop in `run()`, the cadence in `on_state()`, error-correction handshake, concrete subclass tour, pitfalls.
- **`logging-and-progress.md`** — `Logging` / `SubLogging` / `Progress`, the parent walk in `update()`, persistence, common patterns, pitfalls.
- **`inputnum-parsing.md`** — `InputNum::parse` classification; the `(<num>)/F` cofactor-divisor form; `Phi`/`Quad`/`Hex` and auto-detected algebraic factoring; the negative-`Giant` factorial/primorial encoding; `mod`/fingerprinting; `is_half_factored`.
- **`state-serialization.md`** — the `Writer`/`Reader`/`File` on-disk format; the `TaskState` record and the checkpoint lifecycle; `Giant`/`SerializedGWNum` encoding; atomic write + `.md5` sidecar; `unique_fingerprint`. Each consumer owns its own `TYPE` registry (PRST's is documented with PRST, in `checkpoints.md`).
- **`arithmetic-foundation.md`** — the bignum substrate: the `Giant` heap model and its `Giants`/`GMP`/`GW` backends; `GWState::setup`; `GWArithmetic`/`GWNum`; `SerializedGWNum`; the `ReliableGWArithmetic` round-off escalation; and the `integer.*` number-theory helpers + prime sieve.
- **`curves-and-polynomials.md`** — the specialized arithmetic beside the bignum layer: full vs. differential group split; the five concrete arithmetics; `get_NAF_W`/`precomputed_DAC_S_d` chain construction; the liveness map of which paths a consumer (PRST) actually uses.
- **`config-dsl.md`** — the `Config` option-parser DSL: the runtime `ConfigObject` tree vs. the CRTP `*Setup` builder; the `delim` matching convention; the builder verbs; the `parse_args` walk + group fixpoint + exclusive first-match; `.ini` parsing.
- **`container-format.md`** — the `.pack` bundle format: the newline-delimited-JSON-record on-disk layout; `Packer`/`FileContainer`/`FilePacked`; codecs, MD5, and the `container_error` recovery codes.

### Consumer (application) deep-dives — in the `patnashev/prst` repo (`docs/`)

How PRST drives this library is documented with the application: `run-hierarchy.md`, `proof-system.md`,
`abc-batch.md`, `boinc-and-net.md`, `test-harness.md`, `exponentiation-algorithms.md`,
`math-and-theorems.md`. Start from that repo's `lay-of-the-land.md`.

## Coverage status

Every in-scope framework subsystem has a deep-dive (the eight above). Intentionally scoped out: the
prebuilt third-party libraries (`gwnum/` — Woltman's FFT bignum — and `gmp/`), understood only at their
interface; `bow/bow.cpp`, which is BOINC-only build glue relevant only under a consumer's `BOINC` build;
and the fenced low-level details (`group.cpp`'s curve/poly formulas and `container.cpp`'s writer
byte-framing), whose docs cover the interface and on-disk layout and point to the code for the rest.

A caution for maintainers: these docs embed verified line numbers and "this is the only live caller"
claims that reflect the code at the time of writing — re-grep before trusting any such claim. Some claims
about *who consumes* a given primitive depend on the consumer; citations into `prst/src/...` files point
at the consumer's code in the `patnashev/prst` repo, not this one.
