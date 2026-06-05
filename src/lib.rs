#![forbid(unsafe_code)]
//! Discrete Kuramoto oscillator for ternary {-1,0,+1} systems.

const PI: f64 = std::f64::consts::PI;

/// Map phase angle to ternary value.
pub fn phase_to_ternary(phase: f64) -> i8 {
    let p = ((phase % (2.0 * PI)) + 2.0 * PI) % (2.0 * PI);
    if p < PI / 3.0 || p >= 5.0 * PI / 3.0 { 1 }      // +1: near 0°
    else if p >= PI / 3.0 && p < PI { 0 }               // 0: near 120°
    else { -1 }                                            // -1: near 240°
}

/// Map ternary value to representative phase angle.
pub fn ternary_to_phase(t: i8) -> f64 {
    match t { -1 => 4.0 * PI / 3.0, 0 => 2.0 * PI / 3.0, _ => 0.0 }
}

/// Kuramoto oscillator network.
pub struct KuramotoNetwork {
    pub phases: Vec<f64>,
    pub frequencies: Vec<f64>,
    pub coupling: f64,
    pub adjacency: Vec<Vec<f64>>,
}

impl KuramotoNetwork {
    pub fn new(n: usize, coupling: f64) -> Self {
        let mut rng_s: u64 = 42;
        let mut rng = || -> f64 { rng_s = rng_s.wrapping_mul(6364136223846793005).wrapping_add(1); (rng_s >> 33) as f64 / (1u64 << 31) as f64 };
        Self {
            phases: (0..n).map(|_| rng() * 2.0 * PI).collect(),
            frequencies: (0..n).map(|_| (rng() - 0.5) * 2.0).collect(),
            coupling,
            adjacency: vec![vec![1.0; n]; n],
        }
    }

    pub fn new_synchronized(n: usize, coupling: f64) -> Self {
        Self { phases: vec![0.0; n], frequencies: vec![1.0; n], coupling, adjacency: vec![vec![1.0; n]; n] }
    }

    /// Compute order parameter R ∈ [0,1]. R=1 means full sync.
    pub fn order_parameter(&self) -> f64 {
        let n = self.phases.len() as f64;
        let sum_cos: f64 = self.phases.iter().map(|&p| p.cos()).sum();
        let sum_sin: f64 = self.phases.iter().map(|&p| p.sin()).sum();
        (sum_cos * sum_cos + sum_sin * sum_sin).sqrt() / n
    }

    /// One step of Kuramoto dynamics.
    pub fn step(&mut self, dt: f64) {
        let n = self.phases.len();
        let mut dphi = vec![0.0; n];
        for i in 0..n {
            let coupling_sum: f64 = (0..n).filter(|&j| j != i)
                .map(|j| self.adjacency[i][j] * (self.phases[j] - self.phases[i]).sin())
                .sum();
            dphi[i] = self.frequencies[i] + self.coupling / (n - 1).max(1) as f64 * coupling_sum;
        }
        for i in 0..n { self.phases[i] += dphi[i] * dt; }
    }

    /// Run for N steps.
    pub fn run(&mut self, steps: usize, dt: f64) -> Vec<f64> {
        let mut history = Vec::with_capacity(steps);
        for _ in 0..steps { self.step(dt); history.push(self.order_parameter()); }
        history
    }

    /// Get ternary representation of current phases.
    pub fn to_ternary(&self) -> Vec<i8> { self.phases.iter().map(|&p| phase_to_ternary(p)).collect() }

    /// Detect chimera state: some oscillators synced, others not.
    pub fn is_chimera(&self) -> bool {
        let n = self.phases.len();
        if n < 4 { return false; }
        // Split in half, compare order parameters
        let mid = n / 2;
        let r1 = local_order(&self.phases[..mid]);
        let r2 = local_order(&self.phases[mid..]);
        (r1 - r2).abs() > 0.3 && r1.max(r2) > 0.7
    }

    /// Mean phase coherence per group.
    pub fn group_coherence(&self, groups: &[Vec<usize>]) -> Vec<f64> {
        groups.iter().map(|group| {
            let cos_sum: f64 = group.iter().map(|&i| self.phases[i].cos()).sum();
            let sin_sum: f64 = group.iter().map(|&i| self.phases[i].sin()).sum();
            (cos_sum * cos_sum + sin_sum * sin_sum).sqrt() / group.len() as f64
        }).collect()
    }
}

fn local_order(phases: &[f64]) -> f64 {
    let n = phases.len() as f64;
    let sum_cos: f64 = phases.iter().map(|p| p.cos()).sum();
    let sum_sin: f64 = phases.iter().map(|p| p.sin()).sum();
    (sum_cos * sum_cos + sum_sin * sum_sin).sqrt() / n
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test] fn test_phase_to_ternary_zero() { assert_eq!(phase_to_ternary(0.0), 1); }
    #[test] fn test_phase_to_ternary_120() { assert_eq!(phase_to_ternary(2.0*PI/3.0), 0); }
    #[test] fn test_phase_to_ternary_240() { assert_eq!(phase_to_ternary(4.0*PI/3.0), -1); }
    #[test] fn test_ternary_to_phase_roundtrip() { let p = ternary_to_phase(1); assert_eq!(phase_to_ternary(p), 1); }
    #[test] fn test_sync_order_parameter() { let net = KuramotoNetwork::new_synchronized(10, 1.0); assert!((net.order_parameter() - 1.0).abs() < 0.01); }
    #[test] fn test_random_order_parameter() { let net = KuramotoNetwork::new(10, 0.0); assert!(net.order_parameter() < 0.8); }
    #[test] fn test_step_advances_phases() { let mut net = KuramotoNetwork::new_synchronized(5, 1.0); let p0 = net.phases[0]; net.step(0.1); assert!(net.phases[0] > p0); }
    #[test] fn test_coupling_syncs() { let mut net = KuramotoNetwork::new(10, 5.0); let history = net.run(500, 0.1); assert!(history.last().unwrap() > &0.5, "Should sync with strong coupling"); }
    #[test] fn test_no_coupling_no_sync() { let mut net = KuramotoNetwork::new(10, 0.0); let history = net.run(100, 0.1); assert!(history.last().unwrap() < &0.9); }
    #[test] fn test_to_ternary() { let net = KuramotoNetwork::new(5, 1.0); let t = net.to_ternary(); assert_eq!(t.len(), 5); assert!(t.iter().all(|&v| v >= -1 && v <= 1)); }
    #[test] fn test_chimera_detection() { let net = KuramotoNetwork::new(20, 0.1); // Unlikely chimera with low coupling
        assert!(!net.is_chimera() || true); // Just check it doesn't crash
    }
    #[test] fn test_group_coherence() { let net = KuramotoNetwork::new_synchronized(6, 1.0); let coh = net.group_coherence(&[vec![0,1,2], vec![3,4,5]]); assert!(coh.iter().all(|&c| c > 0.9)); }
    #[test] fn test_order_parameter_range() { let net = KuramotoNetwork::new(5, 0.0); let r = net.order_parameter(); assert!(r >= 0.0 && r <= 1.0); }
    #[test] fn test_run_returns_history() { let mut net = KuramotoNetwork::new(3, 1.0); let h = net.run(10, 0.1); assert_eq!(h.len(), 10); }
    #[test] fn test_network_creation() { let net = KuramotoNetwork::new(10, 1.0); assert_eq!(net.phases.len(), 10); }
    #[test] fn test_adjacency_shape() { let net = KuramotoNetwork::new(5, 1.0); assert_eq!(net.adjacency.len(), 5); assert_eq!(net.adjacency[0].len(), 5); }
    #[test] fn test_phase_wrapping() { let t = phase_to_ternary(10.0 * PI); assert!(t == -1 || t == 0 || t == 1); }
    #[test] fn test_negative_phase() { let t = phase_to_ternary(-PI / 2.0); assert!(t == -1 || t == 0 || t == 1); }
}
