"""
Enhanced DNA Sequence Assembler
Robust assembly with reverse complement support and error handling
"""

import numpy as np
from Bio.Align import PairwiseAligner
from collections import defaultdict
from tkinter import messagebox
import sys

class EnhancedAssembly:
    def __init__(self, app):
        self.app = app
        self.similarity_threshold = 0.65  # Lowered for better sensitivity
        self.overlap_threshold = 20       # More lenient overlap requirement
        self.quality_threshold = 20       # Minimum average quality score
        self.min_contig_length = 100      # Minimum contig length to report
        self.max_alignments = 1000        # Prevent alignment explosion
        self.aligner = self._configure_aligner()

    def _configure_aligner(self):
        """Set up optimal alignment parameters"""
        aligner = PairwiseAligner()
        aligner.mode = "local"            # Crucial for overlap detection
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -1.5     # Stricter gap penalties
        aligner.extend_gap_score = -0.5
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0
        return aligner

    def _reverse_complement(self, seq):
        """Generate reverse complement of DNA sequence"""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                     'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
        return ''.join([complement.get(base, 'N') for base in seq[::-1]])

    def _calculate_similarity(self, s1, s2):
        """Calculate sequence similarity score"""
        if not s1 or not s2:
            return 0.0
            
        try:
            alignments = self.aligner.align(s1, s2)[:self.max_alignments]
            if not alignments:
                return 0.0
                
            align = alignments[0]
            aligned1, aligned2 = self._convert_aligned_to_strings(s1, s2, align)
            matches = sum(1 for a, b in zip(aligned1, aligned2) if a == b and a != "-")
            return matches / max(len(s1), len(s2))
        except:
            return 0.0

    def assemble_sequences(self, sequence_names):
        """Main assembly workflow with comprehensive error handling"""
        try:
            # Initial validation
            if len(sequence_names) < 2:
                self._log("Insufficient sequences for assembly")
                return {}

            # Prepare sequences
            sequences = self._prepare_sequences(sequence_names)
            sequences = self._validate_sequences(sequences)
            
            if len(sequences) < 2:
                self._log("Not enough valid sequences after filtering")
                return {}

            self._update_progress(10, "Analyzing sequences...")
            
            # Core assembly steps
            similarity_matrix = self._calculate_similarity_matrix(sequences)
            overlap_data = self._find_overlaps(sequences)
            clusters = self._cluster_sequences(sequences, similarity_matrix, overlap_data)
            
            # Primary assembly attempt
            contigs = self._assemble_clusters(clusters, sequences)
            
            # Fallback assembly if no contigs
            if not contigs:
                self._log("No contigs from clustering, trying fallback assembly")
                contigs = self._fallback_assembly(sequences)
            
            self._update_progress(100, "Assembly complete!")
            return self._filter_contigs(contigs)

        except Exception as e:
            self._show_error(f"Assembly failed: {str(e)}")
            return {}

    def _prepare_sequences(self, sequence_names):
        """Load and preprocess sequences"""
        sequences = {}
        for name in sequence_names:
            if name not in self.app.sequences:
                continue
                
            data = self.app.sequences[name]
            seq = data.get("edited_seq") or str(data["seq_record"].seq)
            is_reverse = data.get("is_reverse", False)
            
            sequences[name] = {
                "processed_seq": self._reverse_complement(seq) if is_reverse else seq,
                "processed_qual": data["quality_scores"][::-1] if is_reverse and isinstance(data.get("quality_scores"), np.ndarray) else data.get("quality_scores"),
                "plot_data": data.get("plot_data", {}),
                "is_reverse": is_reverse,
                "original_length": len(seq)
            }
        return sequences

    def _validate_sequences(self, sequences):
        """Filter out low-quality sequences"""
        valid = {}
        for name, data in sequences.items():
            seq = data["processed_seq"]
            qual = data["processed_qual"]
            
            # Length check
            if len(seq) < self.overlap_threshold:
                continue
                
            # Quality check
            if qual is not None and np.mean(qual) < self.quality_threshold:
                continue
                
            # Complexity check (simple N-count)
            if seq.count('N') / len(seq) > 0.2:
                continue
                
            valid[name] = data
        return valid

    def _calculate_similarity_matrix(self, sequences):
        """Compute pairwise sequence similarities"""
        names = list(sequences.keys())
        n = len(names)
        matrix = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i+1, n):
                s1 = sequences[names[i]]["processed_seq"]
                s2 = sequences[names[j]]["processed_seq"]
                matrix[i][j] = matrix[j][i] = self._calculate_similarity(s1, s2)
        
        return matrix

    def _find_overlaps(self, sequences):
        """Find all significant overlaps between sequences"""
        names = list(sequences.keys())
        overlaps = {}
        
        for i, n1 in enumerate(names):
            for j, n2 in enumerate(names):
                if i >= j:
                    continue
                    
                s1 = sequences[n1]["processed_seq"]
                s2 = sequences[n2]["processed_seq"]
                overlap = self._find_best_overlap(s1, s2)
                
                if overlap["length"] >= self.overlap_threshold:
                    overlaps[(n1, n2)] = overlap
                    
        return overlaps

    def _find_best_overlap(self, s1, s2):
        """Find the best overlap between two sequences"""
        fwd = self._find_overlap(s1, s2)
        rev = self._find_overlap(s2, s1)
        return max(fwd, rev, key=lambda x: (x["length"], x["similarity"]))

    def _find_overlap(self, s1, s2):
        """Detect optimal overlap between sequences"""
        best = {"length": 0, "similarity": 0.0, "score": 0}
        min_len = min(len(s1), len(s2))
        
        # Test multiple overlap windows
        for window_start in range(15, min_len, 25):
            window_end = min(window_start + 50, min_len)
            for i in range(window_start, window_end + 1):
                suffix = s1[-i:]
                prefix = s2[:i]
                
                # Quick similarity estimate
                matches = sum(a == b for a, b in zip(suffix, prefix))
                sim = matches / i
                
                if sim >= 0.7:  # Only do full alignment on promising candidates
                    score = self.aligner.score(suffix, prefix)
                    if score > best["score"]:
                        best = {
                            "length": i,
                            "similarity": sim,
                            "score": score
                        }
        return best

    def _cluster_sequences(self, sequences, sim_matrix, overlap_data):
        """Cluster sequences based on similarity and overlaps"""
        names = list(sequences.keys())
        adj = defaultdict(set)
        
        for i in range(len(names)):
            for j in range(i+1, len(names)):
                if (sim_matrix[i][j] >= self.similarity_threshold or
                    (names[i], names[j]) in overlap_data):
                    adj[names[i]].add(names[j])
                    adj[names[j]].add(names[i])
                    
        return self._find_connected_components(adj)

    def _find_connected_components(self, adj):
        """Find all connected components in the similarity graph"""
        clusters = []
        visited = set()
        
        for node in adj:
            if node not in visited:
                cluster = self._bfs(node, adj, visited)
                if len(cluster) >= 2:
                    clusters.append(cluster)
                    
        return clusters

    def _bfs(self, start, adj, visited):
        """Breadth-first search for connected components"""
        queue = [start]
        cluster = []
        
        while queue:
            node = queue.pop(0)
            if node not in visited:
                visited.add(node)
                cluster.append(node)
                queue.extend(adj[node] - visited)
                
        return cluster

    def _assemble_clusters(self, clusters, sequences):
        """Assemble each cluster into contigs"""
        contigs = {}
        
        for i, cluster in enumerate(clusters):
            contig_name = f"Contig_{i+1}"
            contig_data = self._assemble_cluster(cluster, sequences)
            
            if contig_data and len(contig_data["consensus_sequence"]) >= self.min_contig_length:
                contigs[contig_name] = contig_data
                
        return contigs

    def _assemble_cluster(self, cluster, sequences):
        """Assemble a single cluster of sequences"""
        if len(cluster) < 2:
            return None
            
        order = self._determine_assembly_order(cluster, sequences)
        seq, qual, trace = self._initialize_assembly(sequences[order[0]])
        
        for name in order[1:]:
            seq, qual, trace = self._merge_sequences(
                seq, qual, trace,
                sequences[name]["processed_seq"],
                sequences[name]["processed_qual"],
                sequences[name]["plot_data"]
            )
        
        return self._create_contig_record(seq, qual, trace, cluster, order)

    def _determine_assembly_order(self, cluster, sequences):
        """Determine optimal assembly order for a cluster"""
        forward = [name for name in cluster if not sequences[name]["is_reverse"]]
        reverse = [name for name in cluster if sequences[name]["is_reverse"]]
        
        def sort_key(name):
            seq_len = len(sequences[name]["processed_seq"])
            qual = sequences[name]["processed_qual"]
            qual_score = np.mean(qual) if qual is not None else 0
            return (seq_len, qual_score)  # Prefer longer, higher quality sequences
            
        return sorted(forward, key=sort_key, reverse=True) + \
               sorted(reverse, key=sort_key, reverse=True)

    def _initialize_assembly(self, sequence_data):
        """Initialize assembly with first sequence"""
        return (
            sequence_data["processed_seq"],
            list(sequence_data["processed_qual"]) if sequence_data["processed_qual"] is not None else None,
            {b: list(sequence_data["plot_data"].get(b, [0]*len(sequence_data["processed_seq"]))) for b in "ATCG"}
        )

    def _merge_sequences(self, seq1, qual1, plot1, seq2, qual2, plot2):
        """Merge two sequences using alignment"""
        try:
            alignments = self.aligner.align(seq1, seq2)[:self.max_alignments]
            if not alignments:
                return seq1, qual1, plot1
                
            align = alignments[0]
            a1, a2 = self._convert_aligned_to_strings(seq1, seq2, align)
            return self._build_merged_sequence(a1, a2, qual1, qual2, plot1, plot2)
            
        except Exception as e:
            self._log(f"Merge failed: {str(e)}")
            return seq1, qual1, plot1

    def _convert_aligned_to_strings(self, s1, s2, alignment):
        """Convert alignment object to aligned strings"""
        aligned1 = []
        aligned2 = []
        s1_ptr = 0
        s2_ptr = 0
        
        for (start1, end1), (start2, end2) in zip(*alignment.aligned):
            # Add gaps before aligned block
            while s1_ptr < start1:
                aligned1.append(s1[s1_ptr])
                aligned2.append("-")
                s1_ptr += 1
                
            while s2_ptr < start2:
                aligned1.append("-")
                aligned2.append(s2[s2_ptr])
                s2_ptr += 1
                
            # Add aligned block
            while s1_ptr < end1 and s2_ptr < end2:
                aligned1.append(s1[s1_ptr])
                aligned2.append(s2[s2_ptr])
                s1_ptr += 1
                s2_ptr += 1
                
        # Add remaining sequence
        while s1_ptr < len(s1):
            aligned1.append(s1[s1_ptr])
            aligned2.append("-")
            s1_ptr += 1
            
        while s2_ptr < len(s2):
            aligned1.append("-")
            aligned2.append(s2[s2_ptr])
            s2_ptr += 1
            
        return "".join(aligned1), "".join(aligned2)

    def _build_merged_sequence(self, a1, a2, qual1, qual2, plot1, plot2):
        """Construct merged sequence from alignment"""
        merged_seq = []
        merged_qual = []
        merged_plot = {b: [] for b in "ATCG"}
        
        for b1, b2 in zip(a1, a2):
            if b1 == "-":
                # Take from sequence 2
                merged_seq.append(b2)
                merged_qual.append(self._get_quality(qual2))
                for base in "ATCG":
                    merged_plot[base].append(self._get_plot_value(plot2, base))
                    
            elif b2 == "-":
                # Take from sequence 1
                merged_seq.append(b1)
                merged_qual.append(self._get_quality(qual1))
                for base in "ATCG":
                    merged_plot[base].append(self._get_plot_value(plot1, base))
                    
            else:
                # Merge matching positions
                merged_seq.append(b1 if b1 == b2 else b1)
                q1 = self._get_quality(qual1)
                q2 = self._get_quality(qual2)
                merged_qual.append((q1 + q2) // 2)
                
                for base in "ATCG":
                    v1 = self._get_plot_value(plot1, base)
                    v2 = self._get_plot_value(plot2, base)
                    merged_plot[base].append((v1 + v2) // 2)
                    
        return "".join(merged_seq), merged_qual, merged_plot

    def _get_quality(self, qual):
        """Safely get next quality score"""
        return qual.pop(0) if qual and len(qual) > 0 else 20

    def _get_plot_value(self, plot, base):
        """Safely get next plot value"""
        return plot.get(base, [0]).pop(0) if plot and base in plot and len(plot[base]) > 0 else 0

    def _create_contig_record(self, seq, qual, trace, members, order):
        """Create final contig data structure"""
        return {
            "consensus_sequence": seq,
            "quality_scores": np.array(qual) if qual is not None else None,
            "plot_data": {b: np.array(v) for b, v in trace.items()},
            "members": members,
            "assembly_order": order,
            "metadata": {
                "num_members": len(members),
                "length": len(seq),
                "mean_quality": np.mean(qual) if qual is not None else None
            },
            "file_type": "contig"
        }

    def _fallback_assembly(self, sequences):
        """Attempt assembly when normal clustering fails"""
        contigs = {}
        used = set()
        
        for name1, data1 in sequences.items():
            if name1 in used:
                continue
                
            seq = data1["processed_seq"]
            qual = data1["processed_qual"]
            trace = data1["plot_data"]
            members = [name1]
            
            for name2, data2 in sequences.items():
                if name2 in used or name1 == name2:
                    continue
                    
                overlap = self._find_best_overlap(seq, data2["processed_seq"])
                
                if overlap["length"] >= 15:  # Very lenient threshold
                    seq, qual, trace = self._merge_sequences(
                        seq, qual, trace,
                        data2["processed_seq"],
                        data2["processed_qual"],
                        data2["plot_data"]
                    )
                    members.append(name2)
                    used.add(name2)
            
            if len(seq) >= self.min_contig_length:
                contig_name = f"Contig_{len(contigs)+1}"
                contigs[contig_name] = self._create_contig_record(
                    seq, qual, trace, members, members
                )
        
        return contigs

    def _filter_contigs(self, contigs):
        """Filter contigs by quality/length"""
        return {
            name: data for name, data in contigs.items()
            if len(data["consensus_sequence"]) >= self.min_contig_length
        }

    def _update_progress(self, value, message):
        """Update GUI progress display"""
        if hasattr(self.app, 'root') and hasattr(self.app, 'utils'):
            self.app.root.after(0, self.app.utils.update_progress, value, 100, message)

    def _show_error(self, message):
        """Show error message in GUI"""
        if hasattr(self.app, 'root'):
            self.app.root.after(0, messagebox.showerror, "Assembly Error", message)

    def _log(self, message):
        """Log diagnostic messages"""
        print(f"[Assembly] {message}")
        if hasattr(self.app, 'log_message'):
            self.app.log_message(message)
