import numpy as np

class BaselineGibbsSampler:
    def __init__(self, sequences, motif_width=5, pseudocount=1.0, seed=42):
        self.sequences = sequences
        self.num_seqs = len(sequences)
        self.seq_len = len(sequences[0])
        self.W = motif_width
        self.pseudocount = pseudocount
        self.rng = np.random.default_rng(seed)
        
        # DNA/RNA alphabet mapping
        self.alphabet = ['A', 'C', 'G', 'T']
        self.char_to_idx = {c: i for i, c in enumerate(self.alphabet)}
        
        # Convert string sequences to integer matrices for fast NumPy slicing
        self.seq_matrix = np.zeros((self.num_seqs, self.seq_len), dtype=int)
        for i, seq in enumerate(sequences):
            for j, char in enumerate(seq):
                # Map U to T transparently if present
                base = 'T' if char == 'U' else char
                self.seq_matrix[i, j] = self.char_to_idx.get(base, 0) # Default to A if N
                
        # Calculate empirical background frequencies
        unique, counts = np.unique(self.seq_matrix, return_counts=True)
        self.bg_freqs = np.ones(4) * 0.25 # Fallback uniform
        for u, c in zip(unique, counts):
            self.bg_freqs[u] = c / self.seq_matrix.size
            
        # Initialize random start positions (z)
        self.z = self.rng.integers(0, self.seq_len - self.W + 1, size=self.num_seqs)

    def _get_count_matrix(self, exclude_idx=None):
        """Builds the position count matrix, optionally holding out one sequence."""
        counts = np.zeros((4, self.W))
        for i in range(self.num_seqs):
            if i == exclude_idx:
                continue
            start = self.z[i]
            window = self.seq_matrix[i, start:start+self.W]
            for pos, base in enumerate(window):
                counts[base, pos] += 1
        return counts

    def _get_pwm(self, count_matrix):
        """Converts counts to a Position Weight Matrix with Dirichlet pseudocounts."""
        # Add pseudocounts to avoid log(0)
        counts_pseudo = count_matrix + self.pseudocount
        # Normalize columns to sum to 1
        return counts_pseudo / counts_pseudo.sum(axis=0, keepdims=True)

    def step(self):
        """Performs one complete pass of collapsed Gibbs sampling over all sequences."""
        for i in range(self.num_seqs):
            # 1. Hold out sequence i and compute PWM from the rest
            counts = self._get_count_matrix(exclude_idx=i)
            pwm = self._get_pwm(counts)
            
            # 2. Score all valid candidate windows in sequence i
            valid_starts = self.seq_len - self.W + 1
            log_weights = np.zeros(valid_starts)
            
            for j in range(valid_starts):
                window = self.seq_matrix[i, j:j+self.W]
                
                # Log-likelihood of window under the Motif Model (PWM)
                ll_motif = np.sum(np.log(pwm[window, np.arange(self.W)]))
                
                # Log-likelihood of window under the Background Model
                ll_bg = np.sum(np.log(self.bg_freqs[window]))
                
                # The unnormalized log-weight is the log-likelihood ratio (plus prior)
                # For baseline, prior is uniform, so log(prior) is constant and omitted
                log_weights[j] = ll_motif - ll_bg
            
            # 3. Sample a new start position proportional to the weights
            # Subtract max for numerical stability before exponentiating
            weights = np.exp(log_weights - np.max(log_weights))
            probs = weights / np.sum(weights)
            
            self.z[i] = self.rng.choice(valid_starts, p=probs)

    def get_consensus(self):
        """Returns the consensus string of the current PWM."""
        pwm = self._get_pwm(self._get_count_matrix())
        consensus_idx = np.argmax(pwm, axis=0)
        return "".join([self.alphabet[i] for i in consensus_idx])
    
class StructureAwareGibbsSampler(BaselineGibbsSampler):
    def __init__(self, sequences, unpaired_probs, alpha=1.0, motif_width=5, pseudocount=1.0, seed=42):
        # Initialize the baseline sequence matrices and background frequencies
        super().__init__(sequences, motif_width, pseudocount, seed)
        
        # Store the structural prior matrix and weighting hyperparameter
        self.unpaired_probs = unpaired_probs
        self.alpha = alpha

    def step(self):
        """Executes one round of collapsed Gibbs sampling, incorporating structural priors."""
        for i in range(self.num_seqs):
            counts = self._get_count_matrix(exclude_idx=i)
            pwm = self._get_pwm(counts)
            
            valid_starts = self.seq_len - self.W + 1
            log_weights = np.zeros(valid_starts)
            
            for j in range(valid_starts):
                window = self.seq_matrix[i, j:j+self.W]
                
                # Sequence Likelihoods
                ll_motif = np.sum(np.log(pwm[window, np.arange(self.W)]))
                ll_bg = np.sum(np.log(self.bg_freqs[window]))
                
                # Structural Prior: Proportional to the mean unpaired probability of the window
                mean_unpaired = np.mean(self.unpaired_probs[i, j:j+self.W])
                
                # Add epsilon (1e-9) to prevent log(0) domain errors for perfectly paired windows
                ll_prior = self.alpha * np.log(mean_unpaired + 1e-9)
                
                # Combine sequence likelihood ratio with the structural prior
                log_weights[j] = ll_motif - ll_bg + ll_prior
            
            # Numeric stability subtraction before exponentiation
            weights = np.exp(log_weights - np.max(log_weights))
            probs = weights / np.sum(weights)
            
            self.z[i] = self.rng.choice(valid_starts, p=probs)