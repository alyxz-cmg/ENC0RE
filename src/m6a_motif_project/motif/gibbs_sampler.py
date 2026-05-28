import numpy as np

class BaselineGibbsSampler:
    def __init__(self, sequences, motif_width=5, pseudocount=1.0, seed=42):
        self.sequences = sequences
        self.num_seqs = len(sequences)
        self.seq_len = len(sequences[0])
        self.W = motif_width
        self.pseudocount = pseudocount
        self.rng = np.random.default_rng(seed)
        self.ll_history = []

        # DNA/RNA alphabet mapping
        self.alphabet = ['A', 'C', 'G', 'T']
        self.char_to_idx = {c: i for i, c in enumerate(self.alphabet)}

        # Convert string sequences to integer matrices for fast NumPy slicing
        self.seq_matrix = np.zeros((self.num_seqs, self.seq_len), dtype=int)
        for i, seq in enumerate(sequences):
            for j, char in enumerate(seq):
                base = 'T' if char == 'U' else char
                self.seq_matrix[i, j] = self.char_to_idx.get(base, 0)  # Default A if N

        # Empirical background frequencies
        unique, counts = np.unique(self.seq_matrix, return_counts=True)
        self.bg_freqs = np.ones(4) * 0.25
        for u, c in zip(unique, counts):
            self.bg_freqs[u] = c / self.seq_matrix.size

        # Random initial start positions
        self.z = self.rng.integers(0, self.seq_len - self.W + 1, size=self.num_seqs)

    def _get_count_matrix(self, exclude_idx=None):
        counts = np.zeros((4, self.W))
        for i in range(self.num_seqs):
            if i == exclude_idx:
                continue
            start = self.z[i]
            window = self.seq_matrix[i, start:start + self.W]
            for pos, base in enumerate(window):
                counts[base, pos] += 1
        return counts

    def _get_pwm(self, count_matrix):
        counts_pseudo = count_matrix + self.pseudocount
        return counts_pseudo / counts_pseudo.sum(axis=0, keepdims=True)

    def get_pwm(self):
        """Return the current full-corpus PWM (4 x W)."""
        return self._get_pwm(self._get_count_matrix())

    def step(self):
        for i in range(self.num_seqs):
            counts = self._get_count_matrix(exclude_idx=i)
            pwm = self._get_pwm(counts)

            valid_starts = self.seq_len - self.W + 1
            log_weights = np.zeros(valid_starts)

            for j in range(valid_starts):
                window = self.seq_matrix[i, j:j + self.W]
                ll_motif = np.sum(np.log(pwm[window, np.arange(self.W)]))
                ll_bg = np.sum(np.log(self.bg_freqs[window]))
                log_weights[j] = ll_motif - ll_bg

            weights = np.exp(log_weights - np.max(log_weights))
            probs = weights / np.sum(weights)
            self.z[i] = self.rng.choice(valid_starts, p=probs)

        self.ll_history.append(self.get_log_likelihood())

    def get_consensus(self):
        pwm = self.get_pwm()
        consensus_idx = np.argmax(pwm, axis=0)
        return "".join([self.alphabet[i] for i in consensus_idx])

    def get_log_likelihood(self):
        pwm = self.get_pwm()
        ll_total = 0.0
        for i in range(self.num_seqs):
            start = self.z[i]
            window = self.seq_matrix[i, start:start + self.W]
            ll_total += np.sum(np.log(pwm[window, np.arange(self.W)]))
        return ll_total


class StructureAwareGibbsSampler(BaselineGibbsSampler):
    def __init__(self, sequences, unpaired_probs, alpha=1.0, motif_width=5, pseudocount=1.0, seed=42):
        super().__init__(sequences, motif_width, pseudocount, seed)
        self.unpaired_probs = unpaired_probs
        self.alpha = alpha

    def step(self):
        for i in range(self.num_seqs):
            counts = self._get_count_matrix(exclude_idx=i)
            pwm = self._get_pwm(counts)

            valid_starts = self.seq_len - self.W + 1
            log_weights = np.zeros(valid_starts)

            for j in range(valid_starts):
                window = self.seq_matrix[i, j:j + self.W]
                ll_motif = np.sum(np.log(pwm[window, np.arange(self.W)]))
                ll_bg = np.sum(np.log(self.bg_freqs[window]))

                mean_unpaired = np.mean(self.unpaired_probs[i, j:j + self.W])
                ll_prior = self.alpha * np.log(mean_unpaired + 1e-9)

                log_weights[j] = ll_motif - ll_bg + ll_prior

            weights = np.exp(log_weights - np.max(log_weights))
            probs = weights / np.sum(weights)
            self.z[i] = self.rng.choice(valid_starts, p=probs)

        self.ll_history.append(self.get_log_likelihood())