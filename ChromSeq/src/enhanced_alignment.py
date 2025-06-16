"""
Enhanced alignment module with improved algorithms and visualization
Provides advanced sequence alignment with quality-aware consensus generation
"""

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
import numpy as np
from tkinter import *
from tkinter import ttk, messagebox
from collections import defaultdict
import threading
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

class EnhancedAlignment:
    def __init__(self, app):
        self.app = app
        self.current_alignment = None
        self.alignment_quality = None
        self.alignment_parameters = {
            'match_score': 2,
            'mismatch_score': -1,
            'open_gap_score': -2,
            'extend_gap_score': -0.5,
            'mode': 'global'
        }
        self.colors = {
            'A': '#00CC00', 'T': '#FF0000', 'C': '#0000FF',
            'G': '#000000', 'gap': '#888888'
        }
        self.base_colors = {
            'A': self.colors['A'], 'T': self.colors['T'],
            'C': self.colors['C'], 'G': self.colors['G'],
            'a': self.colors['A'], 't': self.colors['T'],
            'c': self.colors['C'], 'g': self.colors['G'],
            'N': '#FFA500', '-': self.colors['gap']
        }
    
    def align_sequences(self):
        """Enhanced sequence alignment with quality consideration"""
        # Get forward and reverse sequences
        forward_seqs = [(name, data) for name, data in self.app.sequences.items() 
                       if data.get('role') == 'forward']
        reverse_seqs = [(name, data) for name, data in self.app.sequences.items() 
                       if data.get('role') == 'reverse']
        
        if len(forward_seqs) != 1 or len(reverse_seqs) != 1:
            messagebox.showwarning("Warning", 
                                 "Please assign exactly one forward and one reverse sequence")
            return
        
        forward_name, forward_data = forward_seqs[0]
        reverse_name, reverse_data = reverse_seqs[0]
        
        # Show alignment parameters dialog
        if not self.show_alignment_parameters_dialog():
            return
        
        self.app.utils.update_progress(0, 100, "Preparing alignment...")
        
        # Run alignment in background
        self.app.executor.submit(self._perform_enhanced_alignment, 
                               forward_name, forward_data, reverse_name, reverse_data)
    
    def show_alignment_parameters_dialog(self):
        """Show dialog for alignment parameters"""
        dialog = Toplevel(self.app.root)
        dialog.title("Alignment Parameters")
        dialog.geometry("400x350")
        dialog.transient(self.app.root)
        dialog.grab_set()
        
        # Center dialog
        dialog.geometry("+%d+%d" % (self.app.root.winfo_rootx() + 100, 
                                   self.app.root.winfo_rooty() + 100))
        
        # Parameters frame
        params_frame = ttk.LabelFrame(dialog, text="Scoring Parameters", padding="10")
        params_frame.pack(fill=X, padx=20, pady=10)
        
        # Match score
        match_frame = ttk.Frame(params_frame)
        match_frame.pack(fill=X, pady=2)
        ttk.Label(match_frame, text="Match score:").pack(side=LEFT)
        match_var = DoubleVar(value=self.alignment_parameters['match_score'])
        ttk.Spinbox(match_frame, from_=0, to=10, textvariable=match_var, width=10).pack(side=RIGHT)
        
        # Mismatch score
        mismatch_frame = ttk.Frame(params_frame)
        mismatch_frame.pack(fill=X, pady=2)
        ttk.Label(mismatch_frame, text="Mismatch score:").pack(side=LEFT)
        mismatch_var = DoubleVar(value=self.alignment_parameters['mismatch_score'])
        ttk.Spinbox(mismatch_frame, from_=-10, to=0, textvariable=mismatch_var, width=10).pack(side=RIGHT)
        
        # Gap open score
        gap_open_frame = ttk.Frame(params_frame)
        gap_open_frame.pack(fill=X, pady=2)
        ttk.Label(gap_open_frame, text="Gap open score:").pack(side=LEFT)
        gap_open_var = DoubleVar(value=self.alignment_parameters['open_gap_score'])
        ttk.Spinbox(gap_open_frame, from_=-20, to=0, textvariable=gap_open_var, width=10).pack(side=RIGHT)
        
        # Gap extend score
        gap_extend_frame = ttk.Frame(params_frame)
        gap_extend_frame.pack(fill=X, pady=2)
        ttk.Label(gap_extend_frame, text="Gap extend score:").pack(side=LEFT)
        gap_extend_var = DoubleVar(value=self.alignment_parameters['extend_gap_score'])
        ttk.Spinbox(gap_extend_frame, from_=-10, to=0, textvariable=gap_extend_var, width=10).pack(side=RIGHT)
        
        # Alignment mode
        mode_frame = ttk.LabelFrame(dialog, text="Alignment Mode", padding="10")
        mode_frame.pack(fill=X, padx=20, pady=10)
        
        mode_var = StringVar(value=self.alignment_parameters['mode'])
        ttk.Radiobutton(mode_frame, text="Global (Needleman-Wunsch)", 
                       variable=mode_var, value='global').pack(anchor=W)
        ttk.Radiobutton(mode_frame, text="Local (Smith-Waterman)", 
                       variable=mode_var, value='local').pack(anchor=W)
        
        # Quality options
        quality_frame = ttk.LabelFrame(dialog, text="Quality Options", padding="10")
        quality_frame.pack(fill=X, padx=20, pady=10)
        
        use_quality_var = BooleanVar(value=True)
        ttk.Checkbutton(quality_frame, text="Use quality scores in alignment", 
                       variable=use_quality_var).pack(anchor=W)
        
        quality_weight_frame = ttk.Frame(quality_frame)
        quality_weight_frame.pack(fill=X, pady=5)
        ttk.Label(quality_weight_frame, text="Quality weight:").pack(side=LEFT)
        quality_weight_var = DoubleVar(value=0.1)
        ttk.Spinbox(quality_weight_frame, from_=0, to=1, increment=0.1, 
                   textvariable=quality_weight_var, width=10).pack(side=RIGHT)
        
        # Result variable
        result = {'proceed': False}
        
        def apply_parameters():
            self.alignment_parameters = {
                'match_score': match_var.get(),
                'mismatch_score': mismatch_var.get(),
                'open_gap_score': gap_open_var.get(),
                'extend_gap_score': gap_extend_var.get(),
                'mode': mode_var.get(),
                'use_quality': use_quality_var.get(),
                'quality_weight': quality_weight_var.get()
            }
            result['proceed'] = True
            dialog.destroy()
        
        # Buttons
        btn_frame = ttk.Frame(dialog)
        btn_frame.pack(fill=X, padx=20, pady=10)
        
        ttk.Button(btn_frame, text="Align", command=apply_parameters).pack(side=RIGHT, padx=5)
        ttk.Button(btn_frame, text="Cancel", command=dialog.destroy).pack(side=RIGHT, padx=5)
        
        # Wait for dialog to close
        self.app.root.wait_window(dialog)
        return result['proceed']
    
    def _perform_enhanced_alignment(self, forward_name, forward_data, reverse_name, reverse_data):
        """Perform enhanced alignment with quality consideration"""
        try:
            self.app.root.after(0, self.app.utils.update_progress, 10, 100, "Preparing sequences...")
            
            # Get sequences
            forward_seq = forward_data['edited_seq']
            reverse_seq = reverse_data['edited_seq']
            
            if not forward_seq or not reverse_seq:
                self.app.root.after(0, messagebox.showwarning, "Warning", 
                                  "Empty sequences cannot be aligned")
                return
            
            # Get reverse complement of reverse sequence
            reverse_seq_rc = str(Seq(reverse_seq).reverse_complement())
            
            self.app.root.after(0, self.app.utils.update_progress, 30, 100, "Configuring aligner...")
            
            # Configure aligner
            aligner = PairwiseAligner()
            aligner.mode = self.alignment_parameters['mode']
            aligner.match_score = self.alignment_parameters['match_score']
            aligner.mismatch_score = self.alignment_parameters['mismatch_score']
            aligner.open_gap_score = self.alignment_parameters['open_gap_score']
            aligner.extend_gap_score = self.alignment_parameters['extend_gap_score']
            
            self.app.root.after(0, self.app.utils.update_progress, 50, 100, "Performing alignment...")
            
            # Perform alignment
            alignments = aligner.align(forward_seq, reverse_seq_rc)
            
            if not alignments:
                self.app.root.after(0, messagebox.showerror, "Error", 
                                  "No alignment could be found")
                return
            
            # Get best alignment
            best_alignment = alignments[0]
            
            self.app.root.after(0, self.app.utils.update_progress, 70, 100, "Processing alignment...")
            
            # Store alignment result
            # Store alignment result in the app
            self.app.alignment_result = {
                'forward_aligned': str(best_alignment[0]),
                'reverse_aligned': str(best_alignment[1]),
                'forward_data': forward_data,
                'reverse_data': reverse_data,
                'score': best_alignment.score,
                'forward_name': forward_name,
                'reverse_name': reverse_name
            }
            
            
            # Calculate alignment quality metrics
            self._calculate_alignment_quality()
            
            # Generate consensus sequence
            self._generate_consensus()
            
            self.app.root.after(0, self.app.utils.update_progress, 90, 100, "Updating display...")
            
            # Update UI
            self.app.root.after(0, self._display_alignment_results)
            self.app.root.after(0, self.app.utils.reset_progress)
            
        except Exception as e:
            self.app.root.after(0, messagebox.showerror, "Alignment Error", 
                              f"Alignment failed:\\n{str(e)}")
            self.app.root.after(0, self.app.utils.reset_progress)
    
    def _generate_consensus(self):
        """Generate consensus sequence from alignment"""
        if not self.current_alignment:
            return
            
        forward_aligned = self.current_alignment['forward_aligned']
        reverse_aligned = self.current_alignment['reverse_aligned']
        consensus = []
        
        for f_base, r_base in zip(forward_aligned, reverse_aligned):
            if f_base == r_base:
                consensus.append(f_base)
            elif f_base == '-':
                consensus.append(r_base.lower())
            elif r_base == '-':
                consensus.append(f_base.lower())
            else:
                # For mismatches, choose based on quality if available
                f_quality = self._get_base_quality(f_base, 'forward')
                r_quality = self._get_base_quality(r_base, 'reverse')
                
                if f_quality > r_quality:
                    consensus.append(f_base)
                else:
                    consensus.append(r_base.lower())
        
        self.current_alignment['consensus'] = ''.join(consensus)
        self.app.consensus = self.current_alignment['consensus']
    
    def _get_base_quality(self, base, direction):
        """Get quality score for a base (simplified implementation)"""
        if base == '-':
            return 0
        return 20  # Default quality if not available
    
    def _calculate_alignment_quality(self):
        """Calculate quality metrics for the alignment"""
        if not self.current_alignment:
            return
        
        forward_aligned = self.current_alignment['forward_aligned']
        reverse_aligned = self.current_alignment['reverse_aligned']
        forward_data = self.current_alignment['forward_data']
        reverse_data = self.current_alignment['reverse_data']
        
        # Calculate basic statistics
        total_positions = len(forward_aligned)
        matches = sum(1 for f, r in zip(forward_aligned, reverse_aligned) if f == r and f != '-')
        mismatches = sum(1 for f, r in zip(forward_aligned, reverse_aligned) 
                        if f != r and f != '-' and r != '-')
        gaps = sum(1 for f, r in zip(forward_aligned, reverse_aligned) if f == '-' or r == '-')
        
        # Calculate identity percentage
        aligned_length = total_positions - gaps
        identity = (matches / aligned_length * 100) if aligned_length > 0 else 0
        
        # Calculate quality-weighted metrics if quality scores available
        quality_weighted_identity = None
        if (forward_data.get('quality_scores') is not None and 
            reverse_data.get('quality_scores') is not None):
            quality_weighted_identity = self._calculate_quality_weighted_identity(
                forward_aligned, reverse_aligned, 
                forward_data['quality_scores'], reverse_data['quality_scores']
            )
        
        self.alignment_quality = {
            'total_positions': total_positions,
            'matches': matches,
            'mismatches': mismatches,
            'gaps': gaps,
            'identity': identity,
            'quality_weighted_identity': quality_weighted_identity,
            'score': self.current_alignment['score']
        }
    
    def _calculate_quality_weighted_identity(self, forward_aligned, reverse_aligned, 
                                           forward_quality, reverse_quality):
        """Calculate identity weighted by quality scores"""
        weighted_matches = 0
        total_weight = 0
        
        f_pos = 0
        r_pos = 0
        
        for f_base, r_base in zip(forward_aligned, reverse_aligned):
            if f_base != '-' and r_base != '-':
                # Get quality scores for this position
                f_qual = forward_quality[f_pos] if f_pos < len(forward_quality) else 20
                r_qual = reverse_quality[r_pos] if r_pos < len(reverse_quality) else 20
                
                # Use minimum quality as weight
                weight = min(f_qual, r_qual) / 60.0  # Normalize to 0-1
                total_weight += weight
                
                if f_base == r_base:
                    weighted_matches += weight
            
            if f_base != '-':
                f_pos += 1
            if r_base != '-':
                r_pos += 1
        
        return (weighted_matches / total_weight * 100) if total_weight > 0 else 0

    def _display_alignment_results(self):
        """Display alignment results in the main window with consensus and chromatograms"""
        if not self.current_alignment:
            return

        # Create a new window for alignment visualization
        self.alignment_window = Toplevel(self.app.root)
        self.alignment_window.title("Aligned Sequences")
        self.alignment_window.geometry("1200x800")

        # Create container with scrollbars
        container = Frame(self.alignment_window)
        container.pack(fill=BOTH, expand=True)

        # Create canvas with scrollbars
        canvas = Canvas(container)
        scrollbar_y = Scrollbar(container, orient="vertical", command=canvas.yview)
        scrollbar_x = Scrollbar(container, orient="horizontal", command=canvas.xview)
        
        scrollable_frame = Frame(canvas)
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )

        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar_y.set, xscrollcommand=scrollbar_x.set)

        # Pack everything
        scrollbar_x.pack(side=BOTTOM, fill=X)
        scrollbar_y.pack(side=RIGHT, fill=Y)
        canvas.pack(side=LEFT, fill=BOTH, expand=True)

        # Add content to scrollable_frame
        forward_name = self.current_alignment['forward_name']
        reverse_name = self.current_alignment['reverse_name']
        forward_aligned = self.current_alignment['forward_aligned']
        reverse_aligned = self.current_alignment['reverse_aligned']
        forward_data = self.current_alignment['forward_data']
        reverse_data = self.current_alignment['reverse_data']

        # Display alignment header
        header = Label(scrollable_frame, text=f"Alignment of: {forward_name} vs {reverse_name}\nScore: {self.current_alignment['score']:.2f}")
        header.pack(anchor=W)

        if self.alignment_quality:
            qual_label = Label(scrollable_frame, text=f"Identity: {self.alignment_quality['identity']:.1f}%")
            qual_label.pack(anchor=W)

        # Create frames for each sequence with proper scrolling
        frame_height = 200
        current_position = 0
        display_width = 100

        # Ruler frame
        ruler_frame = Frame(scrollable_frame)
        ruler_frame.pack(fill=X, padx=5, pady=2)
        self._create_ruler(ruler_frame, len(forward_aligned))

        # Forward sequence frame
        forward_frame = LabelFrame(scrollable_frame, text="Forward Sequence", height=frame_height)
        forward_frame.pack(fill=X, padx=5, pady=5)
        self._create_sequence_display(forward_frame, forward_data, forward_aligned, "forward")

        # Reverse sequence frame
        reverse_frame = LabelFrame(scrollable_frame, text="Reverse Sequence", height=frame_height)
        reverse_frame.pack(fill=X, padx=5, pady=5)
        self._create_sequence_display(reverse_frame, reverse_data, reverse_aligned, "reverse")

        # Consensus sequence frame if available
        if 'consensus' in self.current_alignment:
            consensus_frame = LabelFrame(scrollable_frame, text="Consensus Sequence", height=frame_height)
            consensus_frame.pack(fill=X, padx=5, pady=5)
            self._create_consensus_display(consensus_frame, self.current_alignment['consensus'])

        # Navigation controls
        nav_frame = Frame(scrollable_frame)
        nav_frame.pack(fill=X, padx=5, pady=5)

        Button(nav_frame, text="<<", command=lambda: self._navigate_alignment(-display_width)).pack(side=LEFT)
        Button(nav_frame, text="<", command=lambda: self._navigate_alignment(-10)).pack(side=LEFT)
        Button(nav_frame, text=">", command=lambda: self._navigate_alignment(10)).pack(side=LEFT)
        Button(nav_frame, text=">>", command=lambda: self._navigate_alignment(display_width)).pack(side=LEFT)

        messagebox.showinfo("Alignment Complete",
                          f"Alignment completed successfully\n"
                          f"Identity: {self.alignment_quality['identity']:.1f}%\n"
                          f"Score: {self.current_alignment['score']:.2f}")

    def _create_ruler(self, parent, length):
        """Create a ruler/position indicator"""
        self.ruler_canvas = Canvas(parent, height=20, bg='white')
        self.ruler_canvas.pack(fill=X)
        
        # Draw the ruler
        self.ruler_canvas.create_line(0, 15, length*10, 15, fill='black')
        
        # Add position markers
        for i in range(0, length, 10):
            if i % 50 == 0:
                # Major tick
                self.ruler_canvas.create_line(i*10, 5, i*10, 15, fill='black')
                self.ruler_canvas.create_text(i*10, 2, text=str(i), anchor=N)
            else:
                # Minor tick
                self.ruler_canvas.create_line(i*10, 10, i*10, 15, fill='black')

    def _create_sequence_display(self, parent, seq_data, aligned_seq, seq_type):
        """Create a display for one sequence type (forward/reverse)"""
        # Create a frame to hold both chromatogram and sequence
        frame = Frame(parent)
        frame.pack(fill=X, expand=True)
        
        # Text widget for sequence display
        text = Text(frame, wrap=NONE, height=10, font=("Courier New", 10))
        scroll_x = Scrollbar(frame, orient=HORIZONTAL, command=text.xview)
        scroll_y = Scrollbar(frame, orient=VERTICAL, command=text.yview)
        text.configure(xscrollcommand=scroll_x.set, yscrollcommand=scroll_y.set)
        
        # Add sequence to text widget
        text.insert(END, aligned_seq)
        text.config(state=DISABLED)
        
        # Pack everything
        scroll_x.pack(side=BOTTOM, fill=X)
        scroll_y.pack(side=RIGHT, fill=Y)
        text.pack(side=LEFT, fill=BOTH, expand=True)       
        # Create a frame for the chromatogram
        fig = Figure(figsize=(10, 2), dpi=100)
        ax = fig.add_subplot(111)
        
        # Plot the chromatogram
        trace_data = seq_data.get('plot_data', {})
        for base in ['A', 'C', 'G', 'T']:
            if base in trace_data:
                ax.plot(trace_data[base], color=self.colors[base], label=base, linewidth=1)
        
        ax.set_ylabel("Intensity")
        ax.legend()
        ax.grid(True, linestyle=':', alpha=0.5)
        
        # Create canvas
        canvas = FigureCanvasTkAgg(fig, master=frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=X, expand=True)
        
        # Create sequence display
        seq_frame = Frame(frame)
        seq_frame.pack(fill=X)
        
        seq_text = Text(
            seq_frame,
            height=1,
            wrap=NONE,
            font=('Courier New', 12),
            padx=5,
            pady=2
        )
        
        # Configure tags for colored bases
        for base, color in self.base_colors.items():
            seq_text.tag_config(base, foreground=color)
        
        # Add aligned sequence
        for base in aligned_seq:
            seq_text.insert(END, base, self.base_colors.get(base, ''))
        
        # Add scrollbar
        scrollbar = Scrollbar(seq_frame, orient="horizontal", command=seq_text.xview)
        seq_text.config(xscrollcommand=scrollbar.set)
        
        seq_text.pack(side=TOP, fill=X)
        scrollbar.pack(side=BOTTOM, fill=X)
        seq_text.config(state=DISABLED)

    def _create_consensus_display(self, parent, consensus_seq):
        """Create a display for the consensus sequence"""
        frame = Frame(parent)
        frame.pack(fill=X, expand=True)
        
        # Create sequence display
        seq_frame = Frame(frame)
        seq_frame.pack(fill=X)
        
        seq_text = Text(
            seq_frame,
            height=1,
            wrap=NONE,
            font=('Courier New', 12),
            padx=5,
            pady=2
        )
        
        # Configure tags for colored bases
        for base, color in self.base_colors.items():
            seq_text.tag_config(base, foreground=color)
        
        # Add consensus sequence
        for base in consensus_seq:
            seq_text.insert(END, base, self.base_colors.get(base, ''))
        
        # Add scrollbar
        scrollbar = Scrollbar(seq_frame, orient="horizontal", command=seq_text.xview)
        seq_text.config(xscrollcommand=scrollbar.set)
        
        seq_text.pack(side=TOP, fill=X)
        scrollbar.pack(side=BOTTOM, fill=X)
        seq_text.config(state=DISABLED)

    def _navigate_alignment(self, delta):
        """Navigate through the alignment (placeholder)"""
        pass

    def get_consensus_sequence(self):
        """Return consensus sequence from last alignment"""
        if not self.current_alignment or 'consensus' not in self.current_alignment:
            return ""
        return self.current_alignment['consensus']

    def get_alignment_statistics(self):
        """Get detailed alignment statistics"""
        if not self.current_alignment or not self.alignment_quality:
            return None
        
        stats = self.alignment_quality.copy()
        
        # Add additional statistics
        forward_aligned = self.current_alignment['forward_aligned']
        reverse_aligned = self.current_alignment['reverse_aligned']
        
        # Calculate gap statistics
        forward_gaps = forward_aligned.count('-')
        reverse_gaps = reverse_aligned.count('-')
        
        # Calculate coverage
        forward_coverage = (len(forward_aligned) - forward_gaps) / len(forward_aligned) * 100
        reverse_coverage = (len(reverse_aligned) - reverse_gaps) / len(reverse_aligned) * 100
        
        stats.update({
            'forward_gaps': forward_gaps,
            'reverse_gaps': reverse_gaps,
            'forward_coverage': forward_coverage,
            'reverse_coverage': reverse_coverage,
            'alignment_length': len(forward_aligned)
        })
        
        return stats
    
    def export_alignment(self):
        """Export alignment to file"""
        if not self.current_alignment:
            messagebox.showwarning("Warning", "No alignment to export")
            return
        
        from tkinter import filedialog
        
        filepath = filedialog.asksaveasfilename(
            title="Export Alignment",
            defaultextension=".aln",
            filetypes=[
                ("Alignment Files", "*.aln"),
                ("Text Files", "*.txt"),
                ("All Files", "*.*")
            ]
        )
        
        if not filepath:
            return
        
        try:
            with open(filepath, 'w') as f:
                f.write("Enhanced Chromas Clone - Alignment Export\n")
                f.write("="*50 + "\n\n")
                
                # Write header information
                f.write(f"Forward sequence: {self.current_alignment['forward_name']}\n")
                f.write(f"Reverse sequence: {self.current_alignment['reverse_name']}\n")
                f.write(f"Alignment score: {self.current_alignment['score']:.2f}\n\n")
                
                if self.alignment_quality:
                    f.write("Alignment Statistics:\n")
                    f.write(f"  Identity: {self.alignment_quality['identity']:.1f}%\n")
                    f.write(f"  Matches: {self.alignment_quality['matches']}\n")
                    f.write(f"  Mismatches: {self.alignment_quality['mismatches']}\n")
                    f.write(f"  Gaps: {self.alignment_quality['gaps']}\n")
                    
                    if self.alignment_quality['quality_weighted_identity'] is not None:
                        f.write(f"  Quality-weighted identity: {self.alignment_quality['quality_weighted_identity']:.1f}%\n")
                
                f.write("\n" + "="*50 + "\n\n")
                
                # Write aligned sequences
                forward_aligned = self.current_alignment['forward_aligned']
                reverse_aligned = self.current_alignment['reverse_aligned']
                
                chunk_size = 60
                for i in range(0, len(forward_aligned), chunk_size):
                    chunk_end = min(i + chunk_size, len(forward_aligned))
                    
                    f.write(f"Position {i+1:>6}\n")
                    f.write(f"Forward:  {forward_aligned[i:chunk_end]}\n")
                    
                    # Alignment indicators
                    indicators = ""
                    for f, r in zip(forward_aligned[i:chunk_end], reverse_aligned[i:chunk_end]):
                        if f == r:
                            indicators += "|"
                        elif f == '-' or r == '-':
                            indicators += " "
                        else:
                            indicators += "."
                    f.write(f"          {indicators}\n")
                    
                    f.write(f"Reverse:  {reverse_aligned[i:chunk_end]}\n\n")
            
            messagebox.showinfo("Success", f"Alignment exported to {filepath}")
            
        except Exception as e:
            messagebox.showerror("Error", f"Could not export alignment:\n{str(e)}")
    
    def clear_alignment(self):
        """Clear current alignment"""
        self.current_alignment = None
        self.alignment_quality = None
        self.app.alignment_result = None
        
        # Clear display
        self.app.sequence_display.config(state=NORMAL)
        self.app.alignment_result.delete(1.0, END)
        self.app.sequence_display.config(state=DISABLED)
        
        # Clear consensus
        self.app.consensus = None
        self.app.consensus_text.config(state=NORMAL)
        self.app.consensus_text.delete(1.0, END)
        self.app.consensus_text.config(state=DISABLED)
