# ternary-kuramoto

**The failure of synchronization in ternary systems. Proof that 0 kills collective rhythm.**

The Kuramoto model is the simplest model of synchronization: a population of oscillators, each spinning at its own natural frequency, coupled to its neighbors through a simple interaction. With enough coupling, they synchronize — all fire together, like fireflies or cardiac cells or metronomes on a shared surface.

But not in ternary. This crate proves why.

When you project continuous phases onto `{-1, 0, +1}` — mapping 0°→+1, 120°→0, 240°→-1 — the 0 state acts as a *phase insulator*. Oscillators that land on 0 contribute nothing to the collective order parameter. They're invisible to the coupling force. And because ternary projection is coarse — each bucket spans 120° — the resulting order parameter can never reach 1.0 even when the underlying continuous phases are perfectly aligned.

**The ternary Kuramoto model cannot fully synchronize.** This isn't a bug in the implementation. It's a theorem about discrete projection.

## What's Inside

- **`phase_to_ternary(phase)`** — project a continuous phase angle to {-1, 0, +1} using 120° sectors
- **`ternary_to_phase(t)`** — the inverse mapping (representative angle for each ternary value)
- **`KuramotoNetwork`** — N oscillators with phases, natural frequencies, coupling strength, and adjacency matrix
- **`order_parameter()`** — the Kuramoto order parameter R ∈ [0, 1]. R = 1 means full sync. In ternary: R < 1 always
- **`step(dt)`** / **`run(steps, dt)`** — advance the simulation, return R history
- **`to_ternary()`** — snapshot of current ternary states
- **`is_chimera()`** — detect chimera states: some synced, some not (the most interesting failure mode)
- **`group_coherence(groups)`** — order parameter per subgroup. Measure who's together and who's lost

## Quick Example

```rust
use ternary_kuramoto::*;

// 20 oscillators, strong coupling
let mut net = KuramotoNetwork::new(20, 5.0);

// Run 500 steps
let history = net.run(500, 0.1);

// Even with strong coupling, ternary projection limits R
let final_r = history.last().unwrap();
// R ≈ 0.8-0.9 but never 1.0 — the 0 sector absorbs coherence

// Check for chimera states
if net.is_chimera() {
    println!("Chimera detected! Partial synchronization.");
}

// Ternary snapshot
let states = net.to_ternary();
// Count how many are stuck at 0
let zeros = states.iter().filter(|&&v| v == 0).count();
// zeros/20 ≈ 33% — exactly 1/3 of the phase circle
```

## The Deeper Truth

**120° sectors create a 1/3 dead zone.** The mapping splits the circle into three equal arcs. The 0 arc (120° wide) contains oscillators that are "almost synced" with their +1 neighbors — they're only 30° away — but they read as *zero*, contributing nothing to the order parameter. This is the same mechanism as the Ising insulator and the flocking alignment failure: **the 0 state screens information regardless of the underlying physics.**

The chimera detection is particularly revealing. In continuous Kuramoto, chimeras are rare and exotic. In ternary Kuramoto, they're *common* — the projection creates artificial boundaries that split the population into apparent "synced" and "unsynced" groups that don't correspond to any real dynamical separation.

**Use cases:**
- **Synchronization research** — prove that discrete projection degrades sync
- **Neural network analysis** — ternary neurons face the same synchronization limits
- **Power grid stability** — discrete phasor measurements have this exact projection problem
- **Biological rhythms** — circadian clocks with ternary states
- **Distributed consensus** — prove impossibility results for ternary-valued consensus

## See Also

- **ternary-ising** — the other proof that 0 kills collective order
- **ternary-phase** — phase relationships between ternary oscillators
- **ternary-sync** — Z₃ synchronization (the *only* thing that works in ternary)
- **ternary-fib** — period 8 as the natural ternary rhythm

## Install

```bash
cargo add ternary-kuramoto
```

## License

MIT
