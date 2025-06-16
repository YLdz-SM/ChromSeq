"""
Enhanced consensus generation with quality-based algorithms
Provides ChromasPro-like consensus functionality with advanced options
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
try:
    import tkinter as tk
    from tkinter import ttk, messagebox
except ImportError:
    import Tkinter as tk  # Python 2
    import ttk
    import tkMessageBox as messagebox
from collections import Counter
import statistics

class EnhancedConsensusGenerator:
    def __init__(self, app):
        self.app = app
        self.consensus_algorithms = {
            'majority': 'Majority Rule',
            'quality_weighted': 'Quality Weighted',
            'threshold': 'Quality Threshold',
            'bayesian': 'Bayesian Consensus',
            'two_step': 'Two-Step Consensus'
        }
        self.current_algorithm = 'quality_weighted'
        self.consensus_parameters = {
            'quality_threshold': 20.0,
            'ambiguity_threshold': 0.3,
            'gap_threshold': 0.5,
            'min_coverage': 2
        }

    def generate_consensus(self):
        """Generate consensus sequence with advanced algorithms"""
        if not self.app.alignment_result:
            messagebox.showwarning("Warning", "Please align sequences first")
            return
        
        # Show consensus parameters dialog
        if not self.show_consensus_parameters_dialog():
            return
        
        try:
            # Get alignment data directly from alignment_result
            forward_aligned = self.app.alignment_result['forward_aligned']
            reverse_aligned = self.app.alignment_result['reverse_aligned']
            forward_data = self.app.alignment_result['forward_data']
            reverse_data = self.app.alignment_result['reverse_data']
            
            # Generate consensus based on selected algorithm
            if self.current_algorithm == 'majority':
                consensus_result = self._majority_consensus(forward_aligned, reverse_aligned)
            elif self.current_algorithm == 'quality_weighted':
                consensus_result = self._quality_weighted_consensus(
                    forward_aligned, reverse_aligned, forward_data, reverse_data)
            elif self.current_algorithm == 'threshold':
                consensus_result = self._threshold_consensus(
                    forward_aligned, reverse_aligned, forward_data, reverse_data)
            elif self.current_algorithm == 'bayesian':
                consensus_result = self._bayesian_consensus(
                    forward_aligned, reverse_aligned, forward_data, reverse_data)
            elif self.current_algorithm == 'two_step':
                consensus_result = self._two_step_consensus(
                    forward_aligned, reverse_aligned, forward_data, reverse_data)
            else:
                raise ValueError(f"Unknown algorithm: {self.current_algorithm}")
            
            # Store and display results
            self.app.consensus = consensus_result['sequence']
            self._display_consensus_results(consensus_result)
            
        except Exception as e:
            messagebox.showerror("Consensus Error", f"Consensus generation failed:\n{str(e)}")
    
    def _two_step_consensus(self, forward_aligned, reverse_aligned, forward_data, reverse_data):
        """Two-step consensus generation:
        1. Generate first consensus by 'Majority Rules'
        2. Generate second consensus using 'Quality Weighted'
        3. Align both consensuses and resolve ambiguities
        """
        # Step 1: Generate majority rules consensus
        majority_consensus = self._majority_consensus(forward_aligned, reverse_aligned)
        
        # Step 2: Generate quality weighted consensus
        quality_consensus = self._quality_weighted_consensus(
            forward_aligned, reverse_aligned, forward_data, reverse_data)
        
        # Step 3: Align and resolve ambiguities
        final_consensus = self._align_and_resolve_ambiguities(
            majority_consensus['sequence'],
            quality_consensus['sequence'],
            forward_aligned, reverse_aligned
        )
        
        # Calculate confidence scores (average of both methods)
        confidence = []
        if majority_consensus.get('confidence') and quality_consensus.get('confidence'):
            for m, q in zip(majority_consensus['confidence'], quality_consensus['confidence']):
                confidence.append((m + q) / 2)
        
        return {
            'sequence': final_consensus,
            'confidence': confidence,
            'algorithm': 'two_step',
            'statistics': {
                'avg_confidence': np.mean(confidence) if confidence else 0,
                'min_confidence': np.min(confidence) if confidence else 0,
                'max_confidence': np.max(confidence) if confidence else 0,
                'majority_avg_conf': np.mean(majority_consensus.get('confidence', [0])),
                'quality_avg_conf': np.mean(quality_consensus.get('confidence', [0]))
            }
        }
    
    def _align_and_resolve_ambiguities(self, seq1, seq2, forward_aligned, reverse_aligned):
        """Align two consensus sequences and resolve ambiguities"""
        final_consensus = []
        iupac_priority = {
            'A': 0, 'C': 1, 'G': 2, 'T': 3,  # Standard bases have highest priority
            'R': 4, 'Y': 5, 'S': 6, 'W': 7, 'K': 8, 'M': 9,  # IUPAC ambiguity codes
            'B': 10, 'D': 11, 'H': 12, 'V': 13,
            'N': 14, '-': 15  # N and gap have lowest priority
        }
        
        for base1, base2 in zip(seq1, seq2):
            if base1 == base2:
                # Perfect match
                final_base = base1
            elif base1 == '-' or base2 == '-':
                # Handle gaps
                final_base = base1 if base2 == '-' else base2
            else:
                # Resolve ambiguity based on priority
                if iupac_priority[base1] < iupac_priority[base2]:
                    final_base = base1
                elif iupac_priority[base2] < iupac_priority[base1]:
                    final_base = base2
                else:
                    # Same priority - use IUPAC code if enabled
                    if self.consensus_parameters.get('use_iupac', True):
                        final_base = self._get_iupac_code([base1, base2])
                    else:
                        final_base = 'N'
            
            final_consensus.append(final_base)
        
        return ''.join(final_consensus)

    def show_consensus_parameters_dialog(self):
        """Show dialog for consensus generation parameters"""
        dialog = tk.Toplevel(self.app.root)
        dialog.title("Consensus Parameters")
        dialog.geometry("450x400")
        dialog.transient(self.app.root)
        dialog.grab_set()
        
        # Center dialog
        dialog.geometry("+%d+%d" % (self.app.root.winfo_rootx() + 100,
                                   self.app.root.winfo_rooty() + 100))
        
        # Algorithm selection
        algo_frame = ttk.LabelFrame(dialog, text="Consensus Algorithm", padding="10")
        algo_frame.pack(fill=tk.X, padx=20, pady=10)
        
        algo_var = tk.StringVar(value=self.current_algorithm)
        for key, name in self.consensus_algorithms.items():
            ttk.Radiobutton(algo_frame, text=name, variable=algo_var, value=key).pack(anchor=tk.W)
        
        # Parameters frame
        params_frame = ttk.LabelFrame(dialog, text="Parameters", padding="10")
        params_frame.pack(fill=tk.X, padx=20, pady=10)
        
        # Quality threshold
        quality_frame = ttk.Frame(params_frame)
        quality_frame.pack(fill=tk.X, pady=2)
        ttk.Label(quality_frame, text="Quality threshold:").pack(side=tk.LEFT)
        quality_var = tk.DoubleVar(value=self.consensus_parameters['quality_threshold'])
        ttk.Spinbox(quality_frame, from_=0, to=60, textvariable=quality_var, width=10).pack(side=tk.RIGHT)
        
        # Ambiguity threshold
        ambig_frame = ttk.Frame(params_frame)
        ambig_frame.pack(fill=tk.X, pady=2)
        ttk.Label(ambig_frame, text="Ambiguity threshold:").pack(side=tk.LEFT)
        ambig_var = tk.DoubleVar(value=self.consensus_parameters['ambiguity_threshold'])
        ttk.Spinbox(ambig_frame, from_=0, to=1, increment=0.1, textvariable=ambig_var, width=10).pack(side=tk.RIGHT)
        
        # Gap threshold
        gap_frame = ttk.Frame(params_frame)
        gap_frame.pack(fill=tk.X, pady=2)
        ttk.Label(gap_frame, text="Gap threshold:").pack(side=tk.LEFT)
        gap_var = tk.DoubleVar(value=self.consensus_parameters['gap_threshold'])
        ttk.Spinbox(gap_frame, from_=0, to=1, increment=0.1, textvariable=gap_var, width=10).pack(side=tk.RIGHT)
        
        # Minimum coverage
        coverage_frame = ttk.Frame(params_frame)
        coverage_frame.pack(fill=tk.X, pady=2)
        ttk.Label(coverage_frame, text="Minimum coverage:").pack(side=tk.LEFT)
        coverage_var = tk.IntVar(value=self.consensus_parameters['min_coverage'])
        ttk.Spinbox(coverage_frame, from_=1, to=10, textvariable=coverage_var, width=10).pack(side=tk.RIGHT)
        
        # Options frame
        options_frame = ttk.LabelFrame(dialog, text="Options", padding="10")
        options_frame.pack(fill=tk.X, padx=20, pady=10)
        
        trim_gaps_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(options_frame, text="Trim terminal gaps", variable=trim_gaps_var).pack(anchor=tk.W)
        
        use_iupac_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(options_frame, text="Use IUPAC ambiguity codes", variable=use_iupac_var).pack(anchor=tk.W)
        
        show_confidence_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(options_frame, text="Show confidence scores", variable=show_confidence_var).pack(anchor=tk.W)
        
        # Result variable
        result = {'proceed': False}
        
        def apply_parameters():
            self.current_algorithm = algo_var.get()
            self.consensus_parameters = {
                'quality_threshold': quality_var.get(),
                'ambiguity_threshold': ambig_var.get(),
                'gap_threshold': gap_var.get(),
                'min_coverage': coverage_var.get(),
                'trim_gaps': trim_gaps_var.get(),
                'use_iupac': use_iupac_var.get(),
                'show_confidence': show_confidence_var.get()
            }
            result['proceed'] = True
            dialog.destroy()
        
        # Buttons
        btn_frame = ttk.Frame(dialog)
        btn_frame.pack(fill=tk.X, padx=20, pady=10)
        
        ttk.Button(btn_frame, text="Generate", command=apply_parameters).pack(side=tk.RIGHT, padx=5)
        ttk.Button(btn_frame, text="Cancel", command=dialog.destroy).pack(side=tk.RIGHT, padx=5)
        
        # Wait for dialog to close
        self.app.root.wait_window(dialog)
        return result['proceed']
    
    def _generate_consensus_background(self):
        """Generate consensus in background thread"""
        try:
            self.app.root.after(0, self.app.utils.update_progress, 20, 100, "Analyzing alignment...")
            
            # Get alignment data
            alignment = self.app.alignment_result
            forward_aligned = alignment['forward_aligned']
            reverse_aligned = alignment['reverse_aligned']
            forward_data = alignment['forward_data']
            reverse_data = alignment['reverse_data']
            
            self.app.root.after(0, self.app.utils.update_progress, 40, 100, "Calculating consensus...")
            
            # Generate consensus based on selected algorithm
            if self.current_algorithm == 'majority':
                consensus_result = self._majority_consensus(forward_aligned, reverse_aligned)
            elif self.current_algorithm == 'quality_weighted':
                consensus_result = self._quality_weighted_consensus(
                    forward_aligned, reverse_aligned, forward_data, reverse_data)
            elif self.current_algorithm == 'threshold':
                consensus_result = self._threshold_consensus(
                    forward_aligned, reverse_aligned, forward_data, reverse_data)
            elif self.current_algorithm == 'bayesian':
                consensus_result = self._bayesian_consensus(
                    forward_aligned, reverse_aligned, forward_data, reverse_data)
            elif self.current_algorithm == 'two_step':
                consensus_result = self._two_step_consensus(
                    forward_aligned, reverse_aligned, forward_data, reverse_data)
            else:
                raise ValueError(f"Unknown algorithm: {self.current_algorithm}")
            
            self.app.root.after(0, self.app.utils.update_progress, 80, 100, "Processing results...")
            
            # Post-process consensus
            if self.consensus_parameters.get('trim_gaps', True):
                consensus_result = self._trim_terminal_gaps(consensus_result)
            
            # Store consensus
            self.app.consensus = consensus_result['sequence']
            self.app.consensus_confidence = consensus_result.get('confidence', [])
            self.app.consensus_statistics = consensus_result.get('statistics', {})
            
            self.app.root.after(0, self.app.utils.update_progress, 100, 100, "Updating display...")
            
            # Update UI
            self.app.root.after(0, self._display_consensus_results, consensus_result)
            self.app.root.after(0, self.app.utils.reset_progress)
            
        except Exception as e:
            self.app.root.after(0, messagebox.showerror, "Consensus Error",
                              f"Consensus generation failed:\n{str(e)}")
            self.app.root.after(0, self.app.utils.reset_progress)
    
    def _majority_consensus(self, forward_aligned, reverse_aligned):
        """Simple majority rule consensus"""
        consensus_seq = []
        confidence_scores = []
        
        for f_base, r_base in zip(forward_aligned, reverse_aligned):
            bases = [b for b in [f_base, r_base] if b != '-']
            
            if not bases:
                consensus_seq.append('-')
                confidence_scores.append(0.0)
            elif len(bases) == 1:
                consensus_seq.append(bases[0])
                confidence_scores.append(0.5)
            else:
                # Count bases
                base_counts = Counter(bases)
                most_common = base_counts.most_common(1)[0]
                consensus_base = most_common[0]
                confidence = most_common[1] / len(bases)
                
                consensus_seq.append(consensus_base)
                confidence_scores.append(confidence)
        
        return {
            'sequence': ''.join(consensus_seq),
            'confidence': confidence_scores,
            'algorithm': 'majority',
            'statistics': {
                'avg_confidence': np.mean(confidence_scores),
                'min_confidence': np.min(confidence_scores),
                'max_confidence': np.max(confidence_scores)
            }
        }
    
    def _quality_weighted_consensus(self, forward_aligned, reverse_aligned,
                                  forward_data, reverse_data):
        """Quality-weighted consensus generation"""
        consensus_seq = []
        confidence_scores = []
        
        # Get quality scores
        forward_quality = forward_data.get('quality_scores')
        reverse_quality = reverse_data.get('quality_scores')
        
        f_pos = 0
        r_pos = 0
        
        for f_base, r_base in zip(forward_aligned, reverse_aligned):
            # Get quality scores for current position
            f_qual = 20.0  # Default quality
            r_qual = 20.0
            
            if f_base != '-' and forward_quality is not None and f_pos < len(forward_quality):
                f_qual = forward_quality[f_pos]
            if r_base != '-' and reverse_quality is not None and r_pos < len(reverse_quality):
                r_qual = reverse_quality[r_pos]
            
            # Determine consensus base
            if f_base == '-' and r_base == '-':
                consensus_base = '-'
                confidence = 0.0
            elif f_base == '-':
                consensus_base = r_base
                confidence = min(r_qual / 60.0, 1.0)
            elif r_base == '-':
                consensus_base = f_base
                confidence = min(f_qual / 60.0, 1.0)
            elif f_base == r_base:
                # Matching bases - high confidence
                consensus_base = f_base
                confidence = min((f_qual + r_qual) / 120.0, 1.0)
            else:
                # Mismatching bases - choose based on quality
                if f_qual > r_qual:
                    consensus_base = f_base
                    confidence = f_qual / (f_qual + r_qual)
                elif r_qual > f_qual:
                    consensus_base = r_base
                    confidence = r_qual / (f_qual + r_qual)
                else:
                    # Equal quality - use ambiguity code if enabled
                    if self.consensus_parameters.get('use_iupac', True):
                        consensus_base = self._get_iupac_code([f_base, r_base])
                    else:
                        consensus_base = 'N'
                    confidence = 0.5
            
            consensus_seq.append(consensus_base)
            confidence_scores.append(confidence)
            
            # Update positions
            if f_base != '-':
                f_pos += 1
            if r_base != '-':
                r_pos += 1
        
        return {
            'sequence': ''.join(consensus_seq),
            'confidence': confidence_scores,
            'algorithm': 'quality_weighted',
            'statistics': {
                'avg_confidence': np.mean(confidence_scores),
                'min_confidence': np.min(confidence_scores),
                'max_confidence': np.max(confidence_scores),
                'high_confidence_positions': sum(1 for c in confidence_scores if c >= 0.8),
                'low_confidence_positions': sum(1 for c in confidence_scores if c < 0.5)
            }
        }
    
    def _threshold_consensus(self, forward_aligned, reverse_aligned,
                           forward_data, reverse_data):
        """Threshold-based consensus generation"""
        quality_threshold = self.consensus_parameters['quality_threshold']
        
        consensus_seq = []
        confidence_scores = []
        
        forward_quality = forward_data.get('quality_scores')
        reverse_quality = reverse_data.get('quality_scores')
        
        f_pos = 0
        r_pos = 0
        
        for f_base, r_base in zip(forward_aligned, reverse_aligned):
            # Get quality scores
            f_qual = 0.0
            r_qual = 0.0
            
            if f_base != '-' and forward_quality is not None and f_pos < len(forward_quality):
                f_qual = forward_quality[f_pos]
            if r_base != '-' and reverse_quality is not None and r_pos < len(reverse_quality):
                r_qual = reverse_quality[r_pos]
            
            # Apply quality threshold
            valid_bases = []
            if f_base != '-' and f_qual >= quality_threshold:
                valid_bases.append((f_base, f_qual))
            if r_base != '-' and r_qual >= quality_threshold:
                valid_bases.append((r_base, r_qual))
            
            if not valid_bases:
                consensus_base = 'N'
                confidence = 0.0
            elif len(valid_bases) == 1:
                consensus_base = valid_bases[0][0]
                confidence = min(valid_bases[0][1] / 60.0, 1.0)
            else:
                # Multiple valid bases
                if valid_bases[0][0] == valid_bases[1][0]:
                    # Same base
                    consensus_base = valid_bases[0][0]
                    confidence = min(sum(qual for _, qual in valid_bases) / 120.0, 1.0)
                else:
                    # Different bases - choose highest quality
                    best_base = max(valid_bases, key=lambda x: x[1])
                    consensus_base = best_base[0]
                    confidence = best_base[1] / sum(qual for _, qual in valid_bases)
            
            consensus_seq.append(consensus_base)
            confidence_scores.append(confidence)
            
            # Update positions
            if f_base != '-':
                f_pos += 1
            if r_base != '-':
                r_pos += 1
        
        return {
            'sequence': ''.join(consensus_seq),
            'confidence': confidence_scores,
            'algorithm': 'threshold',
            'statistics': {
                'avg_confidence': np.mean(confidence_scores),
                'min_confidence': np.min(confidence_scores),
                'max_confidence': np.max(confidence_scores),
                'threshold_used': quality_threshold
            }
        }
    
    def _bayesian_consensus(self, forward_aligned, reverse_aligned,
                          forward_data, reverse_data):
        """Enhanced Bayesian consensus generation with complete probability model"""
        consensus_seq = []
        confidence_scores = []
        
        # Get quality scores and original sequences
        forward_quality = forward_data.get('quality_scores', [])
        reverse_quality = reverse_data.get('quality_scores', [])
        forward_seq = forward_data.get('sequence', '')
        reverse_seq = reverse_data.get('sequence', '')
        
        # Error probability lookup table (Phred scores to probabilities)
        def phred_to_prob(q):
            return 10 ** (-q / 10) if q > 0 else 1.0
        
        f_pos = 0
        r_pos = 0
        
        for f_base, r_base in zip(forward_aligned, reverse_aligned):
            # Initialize base probabilities (A, C, G, T)
            base_probs = np.array([0.0, 0.0, 0.0, 0.0])  # A, C, G, T
            base_map = {'A':0, 'C':1, 'G':2, 'T':3}
            
            # Process forward read
            if f_base != '-':
                # Get original base and quality score
                orig_base = forward_seq[f_pos] if f_pos < len(forward_seq) else 'N'
                f_qual = forward_quality[f_pos] if f_pos < len(forward_quality) else 20.0
                
                # Calculate probabilities
                p_error = phred_to_prob(f_qual)
                p_correct = 1 - p_error
                
                if orig_base in base_map:
                    base_probs[base_map[orig_base]] += np.log(p_correct)
                    # Distribute error probability evenly among other bases
                    for base in ['A', 'C', 'G', 'T']:
                        if base != orig_base:
                            base_probs[base_map[base]] += np.log(p_error/3)
                
                f_pos += 1
            
            # Process reverse read
            if r_base != '-':
                # Get original base and quality score
                orig_base = reverse_seq[r_pos] if r_pos < len(reverse_seq) else 'N'
                r_qual = reverse_quality[r_pos] if r_pos < len(reverse_quality) else 20.0
                
                # Calculate probabilities
                p_error = phred_to_prob(r_qual)
                p_correct = 1 - p_error
                
                if orig_base in base_map:
                    base_probs[base_map[orig_base]] += np.log(p_correct)
                    # Distribute error probability evenly among other bases
                    for base in ['A', 'C', 'G', 'T']:
                        if base != orig_base:
                            base_probs[base_map[base]] += np.log(p_error/3)
                
                r_pos += 1
            
            # Convert from log probabilities back to linear space
            base_probs = np.exp(base_probs - np.max(base_probs))  # Numerical stability
            base_probs = base_probs / np.sum(base_probs)  # Normalize
            
            # Determine consensus base (no IUPAC ambiguity)
            if np.sum(base_probs) > 0:
                consensus_idx = np.argmax(base_probs)
                consensus_base = ['A', 'C', 'G', 'T'][consensus_idx]
                confidence = base_probs[consensus_idx]
            else:
                consensus_base = 'N'
                confidence = 0.0
            
            # Handle gaps - if both reads have gaps at this position
            if f_base == '-' and r_base == '-':
                consensus_base = '-'
                confidence = 0.0
            
            consensus_seq.append(consensus_base)
            confidence_scores.append(confidence)
        
        return {
            'sequence': ''.join(consensus_seq),
            'confidence': confidence_scores,
            'algorithm': 'bayesian',
            'statistics': {
                'avg_confidence': np.mean(confidence_scores),
                'min_confidence': np.min(confidence_scores),
                'max_confidence': np.max(confidence_scores),
                'entropy': np.mean([-np.sum(p*np.log2(p+1e-10)) for p in
                                  [base_probs] if np.sum(base_probs) > 0])
            }
        }
    
    def _get_iupac_code(self, bases):
        """Get IUPAC ambiguity code for a set of bases"""
        base_set = set(bases)
        
        iupac_codes = {
            frozenset(['A', 'G']): 'R',
            frozenset(['C', 'T']): 'Y',
            frozenset(['G', 'C']): 'S',
            frozenset(['A', 'T']): 'W',
            frozenset(['G', 'T']): 'K',
            frozenset(['A', 'C']): 'M',
            frozenset(['C', 'G', 'T']): 'B',
            frozenset(['A', 'G', 'T']): 'D',
            frozenset(['A', 'C', 'T']): 'H',
            frozenset(['A', 'C', 'G']): 'V',
            frozenset(['A', 'C', 'G', 'T']): 'N'
        }
        
        return iupac_codes.get(frozenset(base_set), 'N')
    
    def _trim_terminal_gaps(self, consensus_result):
        """Remove terminal gaps from consensus"""
        sequence = consensus_result['sequence']
        confidence = consensus_result.get('confidence', [])
        
        # Find first and last non-gap positions
        start = 0
        end = len(sequence) - 1
        
        while start < len(sequence) and sequence[start] == '-':
            start += 1
        
        while end >= 0 and sequence[end] == '-':
            end -= 1
        
        if start <= end:
            trimmed_sequence = sequence[start:end+1]
            trimmed_confidence = confidence[start:end+1] if confidence else []
        else:
            trimmed_sequence = ""
            trimmed_confidence = []
        
        result = consensus_result.copy()
        result['sequence'] = trimmed_sequence
        result['confidence'] = trimmed_confidence
        
        return result
    def _display_consensus_results(self, consensus_result):
        """Display consensus results in the UI"""
        # Check if consensus_display exists
        if not hasattr(self.app, 'consensus_display'):
            messagebox.showerror("Display Error", "Consensus display panel not initialized")
            return

        try:
            # Clear previous consensus display
            self.app.consensus_display.config(state=tk.NORMAL)
            self.app.consensus_display.delete(1.0, tk.END)
            
            # Configure text tags for coloring
            self._configure_text_tags()
            
            # Get forward name if available
            forward_name = "Consensus"
            if hasattr(self.app, 'forward_data') and self.app.forward_data:
                forward_name = self.app.forward_data.get('name', 'Consensus')
            
            # Display FASTA header
            fasta_header = f">{forward_name}\n"
            self.app.consensus_display.insert(tk.END, fasta_header)
            
            # Display consensus sequence in FASTA format (60 characters per line)
            sequence = consensus_result['sequence']
            for i in range(0, len(sequence), 60):
                chunk = sequence[i:i+60]
                self.app.consensus_display.insert(tk.END, chunk + "\n")
            
            # Add separator and statistics
            self.app.consensus_display.insert(tk.END, "\n" + "="*50 + "\n\n")
            
            # Display statistics
            statistics = consensus_result.get('statistics', {})
            if statistics:
                stats_text = f"Algorithm: {self.consensus_algorithms[consensus_result['algorithm']]}\n"
                stats_text += f"Average confidence: {statistics.get('avg_confidence', 0):.3f}\n"
                stats_text += f"Length: {len(sequence)} bp\n"
                
                if 'high_confidence_positions' in statistics:
                    stats_text += f"High confidence positions: {statistics['high_confidence_positions']}\n"
                if 'low_confidence_positions' in statistics:
                    stats_text += f"Low confidence positions: {statistics['low_confidence_positions']}\n"
                
                self.app.consensus_display.insert(tk.END, stats_text)
            
            self.app.consensus_display.config(state=tk.DISABLED)
            self.app.consensus_display.see(tk.END)
            
            # Show success message
            messagebox.showinfo("Consensus Generated",
                              f"Consensus sequence generated successfully\n"
                              f"Length: {len(sequence)} bp\n"
                              f"Algorithm: {self.consensus_algorithms[consensus_result['algorithm']]}\n"
                              f"Average confidence: {statistics.get('avg_confidence', 0):.3f}")
        
        except Exception as e:
            messagebox.showerror("Display Error", f"Failed to display consensus: {str(e)}")
            if hasattr(self.app, 'consensus_display'):
                self.app.consensus_display.config(state=tk.DISABLED)

    def _display_consensus_sequence(self, sequence, confidence):
        """Helper method to display the sequence with coloring"""
        chunk_size = 60
        for i in range(0, len(sequence), chunk_size):
            chunk_end = min(i + chunk_size, len(sequence))
            
            # Position line
            self.app.consensus_display.insert(tk.END, f"{i+1:>6}: ")
            
            # Sequence with confidence coloring
            for j in range(i, chunk_end):
                base = sequence[j]
                tag = self._get_confidence_tag(j, confidence, base)
                self.app.consensus_display.insert(tk.END, base, tag)
            
            self.app.consensus_display.insert(tk.END, "\n")
            
            # Confidence line if available
            if confidence and self.consensus_parameters.get('show_confidence', True):
                self.app.consensus_display.insert(tk.END, "       ")
                for j in range(i, chunk_end):
                    if j < len(confidence):
                        conf_char = str(int(confidence[j] * 9))  # 0-9 scale
                    else:
                        conf_char = "?"
                    self.app.consensus_display.insert(tk.END, conf_char, "confidence_scale")
                self.app.consensus_display.insert(tk.END, "\n")
            
            self.app.consensus_display.insert(tk.END, "\n")

    def _configure_text_tags(self):
        """Configure text tags for coloring in the consensus display"""
        # Base colors
        for base, color in self.app.base_colors.items():
            self.app.consensus_display.tag_configure(base, foreground=color)
        
        # Confidence level tags
        self.app.consensus_display.tag_configure("high_confidence", foreground="black")
        self.app.consensus_display.tag_configure("medium_confidence", foreground="orange")
        self.app.consensus_display.tag_configure("low_confidence", foreground="red")
        self.app.consensus_display.tag_configure("confidence_scale", foreground="gray")

    def _get_confidence_tag(self, position, confidence, base):
        """Determine the appropriate tag for a base based on confidence"""
        if position >= len(confidence):
            return base  # Use default base color
        
        conf = confidence[position]
        if conf >= 0.8:
            return "high_confidence"
        elif conf >= 0.5:
            return "medium_confidence"
        else:
            return "low_confidence"

    def copy_consensus(self):
        """Copy consensus sequence to clipboard"""
        if not self.app.consensus:
            messagebox.showwarning("Warning", "No consensus sequence available")
            return
        
        self.app.root.clipboard_clear()
        self.app.root.clipboard_append(self.app.consensus)
        self.app.root.update()
        
        messagebox.showinfo("Copied", "Consensus sequence copied to clipboard")
    
    def export_consensus(self):
        """Export consensus sequence to file"""
        if not self.app.consensus:
            messagebox.showwarning("Warning", "No consensus sequence available")
            return
        
        self.app.file_io.export_fasta(self.app.consensus, "consensus_sequence")
    
    def show_consensus_statistics(self):
        """Show detailed consensus statistics"""
        if not hasattr(self.app, 'consensus_statistics') or not self.app.consensus_statistics:
            messagebox.showwarning("Warning", "No consensus statistics available")
            return
        
        # Create statistics window
        stats_window = tk.Toplevel(self.app.root)
        stats_window.title("Consensus Statistics")
        stats_window.geometry("400x300")
        stats_window.transient(self.app.root)
        
        # Create text widget for statistics
        stats_text = tk.Text(stats_window, wrap=tk.WORD, font=('Consolas', 10))
        scrollbar = ttk.Scrollbar(stats_window, orient=tk.VERTICAL, command=stats_text.yview)
        stats_text.configure(yscrollcommand=scrollbar.set)
        
        stats_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Display statistics
        stats = self.app.consensus_statistics
        stats_content = "Consensus Generation Statistics\n"
        stats_content += "="*40 + "\n\n"
        
        for key, value in stats.items():
            if isinstance(value, float):
                stats_content += f"{key.replace('_', ' ').title()}: {value:.4f}\n"
            else:
                stats_content += f"{key.replace('_', ' ').title()}: {value}\n"
        
        stats_text.insert(1.0, stats_content)
        stats_text.config(state=tk.DISABLED)
