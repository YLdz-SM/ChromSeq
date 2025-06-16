"""
Enhanced sequence handling with interactive editing capabilities
Provides ChromasPro-like editing functionality with undo/redo support
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import os
from tkinter import *
from tkinter import ttk, messagebox, simpledialog
from collections import deque
import copy

class EditOperation:
    """Represents a single edit operation for undo/redo functionality"""
    def __init__(self, operation_type, position, old_value, new_value, sequence_name):
        self.operation_type = operation_type  # 'insert', 'delete', 'substitute'
        self.position = position
        self.old_value = old_value
        self.new_value = new_value
        self.sequence_name = sequence_name
        self.timestamp = None

class EnhancedSequenceHandler:
    def __init__(self, app):
        self.app = app
        self.undo_stack = deque(maxlen=100)  # Limit undo history
        self.redo_stack = deque(maxlen=100)
        self.editing_mode = False
        self.selected_position = None
        self.selected_sequence = None
    
    def toggle_editing_mode(self):
        """Toggle between viewing and editing modes"""
        self.editing_mode = not self.editing_mode
        self.app.edit_mode_var.set(self.editing_mode)
        
        if self.editing_mode:
            self.app.status_bar.config(text="Editing mode enabled - Click on chromatogram to edit bases")
            # Change cursor to indicate editing mode
            self.app.canvas.get_tk_widget().config(cursor="crosshair")
        else:
            self.app.status_bar.config(text="Viewing mode - Editing disabled")
            self.app.canvas.get_tk_widget().config(cursor="")
            self.selected_position = None
            self.selected_sequence = None
        
        # Update visualization to show/hide editing indicators
        self.app.update_visualization()
    
    def handle_chromatogram_click(self, event, sequence_name, position):
        """Handle click events on chromatogram for editing"""
        if not self.editing_mode:
            return
        
        if event.inaxes is None:
            return
        
        # Convert click coordinates to sequence position
        click_position = int(round(event.xdata))
        
        if click_position < 0 or click_position >= len(self.app.sequences[sequence_name]['edited_seq']):
            return
        
        self.selected_position = click_position
        self.selected_sequence = sequence_name
        
        # Show editing options
        self.show_edit_menu(event, sequence_name, click_position)
    
    def show_edit_menu(self, event, sequence_name, position):
        """Show context menu for editing operations"""
        # Create popup menu
        edit_menu = Menu(self.app.root, tearoff=0)
        
        current_base = self.app.sequences[sequence_name]['edited_seq'][position]
        
        # Base substitution options
        for base in ['A', 'T', 'C', 'G', 'N']:
            if base != current_base:
                edit_menu.add_command(
                    label=f"Change to {base}",
                    command=lambda b=base: self.substitute_base(sequence_name, position, b)
                )
        
        edit_menu.add_separator()
        
        # Insertion options
        insert_menu = Menu(edit_menu, tearoff=0)
        for base in ['A', 'T', 'C', 'G', 'N', '-']:
            insert_menu.add_command(
                label=f"Insert {base}",
                command=lambda b=base: self.insert_base(sequence_name, position, b)
            )
        edit_menu.add_cascade(label="Insert Base", menu=insert_menu)
        
        # Deletion option
        edit_menu.add_command(
            label="Delete Base",
            command=lambda: self.delete_base(sequence_name, position)
        )
        
        edit_menu.add_separator()
        
        # Quality-based suggestions
        if self.app.sequences[sequence_name].get('quality_scores') is not None:
            quality_score = self.app.sequences[sequence_name]['quality_scores'][position]
            if quality_score < 20:
                edit_menu.add_command(
                    label=f"Low quality ({quality_score:.1f}) - Suggest N",
                    command=lambda: self.substitute_base(sequence_name, position, 'N')
                )
        
        # Show menu at click position
        try:
            edit_menu.tk_popup(event.x_root, event.y_root)
        finally:
            edit_menu.grab_release()
    
    def substitute_base(self, sequence_name, position, new_base):
        """Substitute a base at the specified position"""
        if sequence_name not in self.app.sequences:
            return
        
        seq_data = self.app.sequences[sequence_name]
        old_sequence = seq_data['edited_seq']
        
        if position < 0 or position >= len(old_sequence):
            return
        
        old_base = old_sequence[position]
        
        # Create new sequence
        new_sequence = old_sequence[:position] + new_base + old_sequence[position+1:]
        
        # Create edit operation for undo
        edit_op = EditOperation('substitute', position, old_base, new_base, sequence_name)
        self.add_to_undo_stack(edit_op)
        
        # Apply the change
        seq_data['edited_seq'] = new_sequence
        
        # Update trace data if needed
        self.update_trace_data_for_edit(sequence_name, 'substitute', position, new_base)
        
        # Update visualization
        self.app.update_visualization()
        
        # Update sequence info
        self.app.utils.update_sequence_info(sequence_name)
        
        self.app.status_bar.config(text=f"Changed {old_base} to {new_base} at position {position+1}")
    
    def insert_base(self, sequence_name, position, new_base):
        """Insert a base at the specified position"""
        if sequence_name not in self.app.sequences:
            return
        
        seq_data = self.app.sequences[sequence_name]
        old_sequence = seq_data['edited_seq']
        
        # Create new sequence
        new_sequence = old_sequence[:position] + new_base + old_sequence[position:]
        
        # Create edit operation for undo
        edit_op = EditOperation('insert', position, '', new_base, sequence_name)
        self.add_to_undo_stack(edit_op)
        
        # Apply the change
        seq_data['edited_seq'] = new_sequence
        
        # Update trace data
        self.update_trace_data_for_edit(sequence_name, 'insert', position, new_base)
        
        # Update visualization
        self.app.update_visualization()
        
        # Update sequence info
        self.app.utils.update_sequence_info(sequence_name)
        
        self.app.status_bar.config(text=f"Inserted {new_base} at position {position+1}")
    
    def delete_base(self, sequence_name, position):
        """Delete a base at the specified position"""
        if sequence_name not in self.app.sequences:
            return
        
        seq_data = self.app.sequences[sequence_name]
        old_sequence = seq_data['edited_seq']
        
        if position < 0 or position >= len(old_sequence):
            return
        
        old_base = old_sequence[position]
        
        # Create new sequence
        new_sequence = old_sequence[:position] + old_sequence[position+1:]
        
        # Create edit operation for undo
        edit_op = EditOperation('delete', position, old_base, '', sequence_name)
        self.add_to_undo_stack(edit_op)
        
        # Apply the change
        seq_data['edited_seq'] = new_sequence
        
        # Update trace data
        self.update_trace_data_for_edit(sequence_name, 'delete', position, '')
        
        # Update visualization
        self.app.update_visualization()
        
        # Update sequence info
        self.app.utils.update_sequence_info(sequence_name)
        
        self.app.status_bar.config(text=f"Deleted {old_base} at position {position+1}")
    
    def update_trace_data_for_edit(self, sequence_name, operation, position, new_base):
        """Update trace data when sequence is edited"""
        seq_data = self.app.sequences[sequence_name]
        
        if 'plot_data' not in seq_data:
            return
        
        trace_data = seq_data['plot_data']
        
        if operation == 'substitute':
            # For substitution, we can try to adjust the trace values
            # This is a simplified approach - in reality, trace editing is complex
            if new_base in trace_data and position < len(trace_data[new_base]):
                # Boost the signal for the new base
                for base in trace_data:
                    if base == new_base:
                        # Increase signal for the new base
                        if position < len(trace_data[base]):
                            trace_data[base][position] = max(trace_data[base][position], 1000)
                    else:
                        # Decrease signal for other bases
                        if position < len(trace_data[base]):
                            trace_data[base][position] = min(trace_data[base][position], 100)
        
        elif operation == 'insert':
            # For insertion, we need to insert a position in all trace arrays
            for base in trace_data:
                if position < len(trace_data[base]):
                    # Insert a new value (interpolated from neighbors)
                    if position > 0 and position < len(trace_data[base]) - 1:
                        new_value = (trace_data[base][position-1] + trace_data[base][position]) // 2
                    else:
                        new_value = trace_data[base][position] if position < len(trace_data[base]) else 0
                    
                    # Convert to list, insert, convert back to array
                    trace_list = trace_data[base].tolist()
                    trace_list.insert(position, new_value)
                    trace_data[base] = np.array(trace_list, dtype=np.int16)
        
        elif operation == 'delete':
            # For deletion, remove the position from all trace arrays
            for base in trace_data:
                if position < len(trace_data[base]):
                    # Convert to list, delete, convert back to array
                    trace_list = trace_data[base].tolist()
                    if position < len(trace_list):
                        del trace_list[position]
                        trace_data[base] = np.array(trace_list, dtype=np.int16)
        
        # Update quality scores and peak positions similarly
        if 'quality_scores' in seq_data and seq_data['quality_scores'] is not None:
            self.update_quality_scores_for_edit(seq_data, operation, position)
        
        if 'peak_positions' in seq_data and seq_data['peak_positions'] is not None:
            self.update_peak_positions_for_edit(seq_data, operation, position)
    
    def update_quality_scores_for_edit(self, seq_data, operation, position):
        """Update quality scores when sequence is edited"""
        quality_scores = seq_data['quality_scores']
        
        if operation == 'insert':
            # Insert a quality score (average of neighbors or default)
            if position > 0 and position < len(quality_scores):
                new_quality = (quality_scores[position-1] + quality_scores[position]) / 2
            else:
                new_quality = 20.0  # Default quality for inserted bases
            
            quality_list = quality_scores.tolist()
            quality_list.insert(position, new_quality)
            seq_data['quality_scores'] = np.array(quality_list, dtype=np.float32)
        
        elif operation == 'delete':
            # Remove quality score at position
            if position < len(quality_scores):
                quality_list = quality_scores.tolist()
                del quality_list[position]
                seq_data['quality_scores'] = np.array(quality_list, dtype=np.float32)
    
    def update_peak_positions_for_edit(self, seq_data, operation, position):
        """Update peak positions when sequence is edited"""
        peak_positions = seq_data['peak_positions']
        
        if operation == 'insert':
            # Insert a peak position (interpolated)
            if position > 0 and position < len(peak_positions):
                new_peak = (peak_positions[position-1] + peak_positions[position]) // 2
            else:
                new_peak = position * 10  # Default spacing
            
            peak_list = peak_positions.tolist()
            peak_list.insert(position, new_peak)
            seq_data['peak_positions'] = np.array(peak_list, dtype=np.int32)
        
        elif operation == 'delete':
            # Remove peak position
            if position < len(peak_positions):
                peak_list = peak_positions.tolist()
                del peak_list[position]
                seq_data['peak_positions'] = np.array(peak_list, dtype=np.int32)
    
    def add_to_undo_stack(self, edit_operation):
        """Add an edit operation to the undo stack"""
        self.undo_stack.append(edit_operation)
        self.redo_stack.clear()  # Clear redo stack when new edit is made
    
    def undo_edit(self):
        """Undo the last edit operation"""
        if not self.undo_stack:
            messagebox.showinfo("Undo", "Nothing to undo")
            return
        
        edit_op = self.undo_stack.pop()
        
        # Reverse the operation
        if edit_op.operation_type == 'substitute':
            self._apply_substitute(edit_op.sequence_name, edit_op.position, edit_op.old_value)
        elif edit_op.operation_type == 'insert':
            self._apply_delete(edit_op.sequence_name, edit_op.position)
        elif edit_op.operation_type == 'delete':
            self._apply_insert(edit_op.sequence_name, edit_op.position, edit_op.old_value)
        
        # Add to redo stack
        self.redo_stack.append(edit_op)
        
        # Update visualization
        self.app.update_visualization()
        self.app.utils.update_sequence_info(edit_op.sequence_name)
        
        self.app.status_bar.config(text=f"Undid {edit_op.operation_type} operation")
    
    def redo_edit(self):
        """Redo the last undone edit operation"""
        if not self.redo_stack:
            messagebox.showinfo("Redo", "Nothing to redo")
            return
        
        edit_op = self.redo_stack.pop()
        
        # Reapply the operation
        if edit_op.operation_type == 'substitute':
            self._apply_substitute(edit_op.sequence_name, edit_op.position, edit_op.new_value)
        elif edit_op.operation_type == 'insert':
            self._apply_insert(edit_op.sequence_name, edit_op.position, edit_op.new_value)
        elif edit_op.operation_type == 'delete':
            self._apply_delete(edit_op.sequence_name, edit_op.position)
        
        # Add back to undo stack
        self.undo_stack.append(edit_op)
        
        # Update visualization
        self.app.update_visualization()
        self.app.utils.update_sequence_info(edit_op.sequence_name)
        
        self.app.status_bar.config(text=f"Redid {edit_op.operation_type} operation")
    
    def _apply_substitute(self, sequence_name, position, new_base):
        """Apply substitution without adding to undo stack"""
        seq_data = self.app.sequences[sequence_name]
        old_sequence = seq_data['edited_seq']
        new_sequence = old_sequence[:position] + new_base + old_sequence[position+1:]
        seq_data['edited_seq'] = new_sequence
        self.update_trace_data_for_edit(sequence_name, 'substitute', position, new_base)
    
    def _apply_insert(self, sequence_name, position, new_base):
        """Apply insertion without adding to undo stack"""
        seq_data = self.app.sequences[sequence_name]
        old_sequence = seq_data['edited_seq']
        new_sequence = old_sequence[:position] + new_base + old_sequence[position:]
        seq_data['edited_seq'] = new_sequence
        self.update_trace_data_for_edit(sequence_name, 'insert', position, new_base)
    
    def _apply_delete(self, sequence_name, position):
        """Apply deletion without adding to undo stack"""
        seq_data = self.app.sequences[sequence_name]
        old_sequence = seq_data['edited_seq']
        new_sequence = old_sequence[:position] + old_sequence[position+1:]
        seq_data['edited_seq'] = new_sequence
        self.update_trace_data_for_edit(sequence_name, 'delete', position, '')
    
    def insert_base_dialog(self):
        """Show dialog for inserting bases at current position"""
        if not self.editing_mode:
            messagebox.showwarning("Warning", "Please enable editing mode first")
            return
        
        if self.selected_position is None or self.selected_sequence is None:
            messagebox.showwarning("Warning", "Please click on a chromatogram position first")
            return
        
        # Create dialog
        dialog = Toplevel(self.app.root)
        dialog.title("Insert Base")
        dialog.geometry("300x150")
        dialog.transient(self.app.root)
        dialog.grab_set()
        
        # Center the dialog
        dialog.geometry("+%d+%d" % (self.app.root.winfo_rootx() + 50, self.app.root.winfo_rooty() + 50))
        
        # Base selection
        ttk.Label(dialog, text="Select base to insert:").pack(pady=10)
        
        base_var = StringVar(value='A')
        base_frame = ttk.Frame(dialog)
        base_frame.pack(pady=5)
        
        for base in ['A', 'T', 'C', 'G', 'N', '-']:
            ttk.Radiobutton(base_frame, text=base, variable=base_var, value=base).pack(side=LEFT, padx=5)
        
        # Position info
        ttk.Label(dialog, text=f"Position: {self.selected_position + 1}").pack(pady=5)
        ttk.Label(dialog, text=f"Sequence: {self.selected_sequence}").pack(pady=5)
        
        # Buttons
        btn_frame = ttk.Frame(dialog)
        btn_frame.pack(pady=10)
        
        def insert_and_close():
            self.insert_base(self.selected_sequence, self.selected_position, base_var.get())
            dialog.destroy()
        
        ttk.Button(btn_frame, text="Insert", command=insert_and_close).pack(side=LEFT, padx=5)
        ttk.Button(btn_frame, text="Cancel", command=dialog.destroy).pack(side=LEFT, padx=5)
    
    def delete_selected_base(self):
        """Delete base at current selected position"""
        if not self.editing_mode:
            messagebox.showwarning("Warning", "Please enable editing mode first")
            return
        
        if self.selected_position is None or self.selected_sequence is None:
            messagebox.showwarning("Warning", "Please click on a chromatogram position first")
            return
        
        # Confirm deletion
        current_base = self.app.sequences[self.selected_sequence]['edited_seq'][self.selected_position]
        if messagebox.askyesno("Confirm Deletion", 
                              f"Delete base '{current_base}' at position {self.selected_position + 1}?"):
            self.delete_base(self.selected_sequence, self.selected_position)
    
    def trim_sequences(self):
        """Trim sequences based on quality scores"""
        selected_items = self.app.file_tree.selection()
        if not selected_items:
            messagebox.showwarning("Warning", "Please select sequences to trim")
            return
        
        # Create trim dialog
        dialog = Toplevel(self.app.root)
        dialog.title("Trim Sequences")
        dialog.geometry("400x300")
        dialog.transient(self.app.root)
        dialog.grab_set()
        
        # Trim parameters
        ttk.Label(dialog, text="Trim Parameters", font=('TkDefaultFont', 12, 'bold')).pack(pady=10)
        
        # Quality threshold
        quality_frame = ttk.Frame(dialog)
        quality_frame.pack(fill=X, padx=20, pady=5)
        ttk.Label(quality_frame, text="Quality threshold:").pack(side=LEFT)
        quality_var = DoubleVar(value=20.0)
        quality_spin = ttk.Spinbox(quality_frame, from_=0, to=60, textvariable=quality_var, width=10)
        quality_spin.pack(side=RIGHT)
        
        # Window size
        window_frame = ttk.Frame(dialog)
        window_frame.pack(fill=X, padx=20, pady=5)
        ttk.Label(window_frame, text="Window size:").pack(side=LEFT)
        window_var = IntVar(value=10)
        window_spin = ttk.Spinbox(window_frame, from_=1, to=50, textvariable=window_var, width=10)
        window_spin.pack(side=RIGHT)
        
        # Trim options
        options_frame = ttk.LabelFrame(dialog, text="Trim Options", padding="10")
        options_frame.pack(fill=X, padx=20, pady=10)
        
        trim_5_var = BooleanVar(value=True)
        trim_3_var = BooleanVar(value=True)
        
        ttk.Checkbutton(options_frame, text="Trim 5' end", variable=trim_5_var).pack(anchor=W)
        ttk.Checkbutton(options_frame, text="Trim 3' end", variable=trim_3_var).pack(anchor=W)
        
        # Preview
        preview_text = Text(dialog, height=6, width=50)
        preview_text.pack(fill=BOTH, expand=True, padx=20, pady=5)
        
        def update_preview():
            preview_text.delete(1.0, END)
            for item in selected_items:
                filename = self.app.file_tree.item(item, 'text')
                if filename in self.app.sequences:
                    seq_data = self.app.sequences[filename]
                    if seq_data.get('quality_scores') is not None:
                        trim_result = self._calculate_trim_positions(
                            seq_data['quality_scores'], 
                            quality_var.get(), 
                            window_var.get(),
                            trim_5_var.get(),
                            trim_3_var.get()
                        )
                        preview_text.insert(END, f"{filename}: {trim_result['start']}-{trim_result['end']} "
                                                f"(length: {trim_result['length']})\\n")
                    else:
                        preview_text.insert(END, f"{filename}: No quality data\\n")
        
        # Update preview when parameters change
        quality_var.trace('w', lambda *args: update_preview())
        window_var.trace('w', lambda *args: update_preview())
        trim_5_var.trace('w', lambda *args: update_preview())
        trim_3_var.trace('w', lambda *args: update_preview())
        
        update_preview()
        
        # Buttons
        btn_frame = ttk.Frame(dialog)
        btn_frame.pack(fill=X, padx=20, pady=10)
        
        def apply_trim():
            for item in selected_items:
                filename = self.app.file_tree.item(item, 'text')
                if filename in self.app.sequences:
                    self._apply_trim_to_sequence(filename, quality_var.get(), window_var.get(),
                                               trim_5_var.get(), trim_3_var.get())
            
            self.app.update_visualization()
            self.app.utils.update_selection_info()
            dialog.destroy()
            messagebox.showinfo("Success", f"Trimmed {len(selected_items)} sequences")
        
        ttk.Button(btn_frame, text="Apply", command=apply_trim).pack(side=RIGHT, padx=5)
        ttk.Button(btn_frame, text="Cancel", command=dialog.destroy).pack(side=RIGHT, padx=5)
    
    def _calculate_trim_positions(self, quality_scores, threshold, window_size, trim_5, trim_3):
        """Calculate trim positions based on quality scores"""
        seq_length = len(quality_scores)
        start_pos = 0
        end_pos = seq_length - 1
        
        if trim_5:
            # Find 5' trim position
            for i in range(seq_length - window_size + 1):
                window_avg = np.mean(quality_scores[i:i+window_size])
                if window_avg >= threshold:
                    start_pos = i
                    break
        
        if trim_3:
            # Find 3' trim position
            for i in range(seq_length - window_size, -1, -1):
                window_avg = np.mean(quality_scores[i:i+window_size])
                if window_avg >= threshold:
                    end_pos = i + window_size - 1
                    break
        
        return {
            'start': start_pos,
            'end': end_pos,
            'length': max(0, end_pos - start_pos + 1)
        }
    
    def _apply_trim_to_sequence(self, filename, threshold, window_size, trim_5, trim_3):
        """Apply trimming to a specific sequence"""
        seq_data = self.app.sequences[filename]
        
        if seq_data.get('quality_scores') is None:
            return
        
        trim_result = self._calculate_trim_positions(
            seq_data['quality_scores'], threshold, window_size, trim_5, trim_3
        )
        
        start_pos = trim_result['start']
        end_pos = trim_result['end']
        
        if start_pos >= end_pos:
            return  # Invalid trim
        
        # Create edit operation for undo
        old_sequence = seq_data['edited_seq']
        new_sequence = old_sequence[start_pos:end_pos+1]
        
        edit_op = EditOperation('trim', 0, old_sequence, new_sequence, filename)
        self.add_to_undo_stack(edit_op)
        
        # Apply trim
        seq_data['edited_seq'] = new_sequence
        
        # Trim associated data
        if 'quality_scores' in seq_data and seq_data['quality_scores'] is not None:
            seq_data['quality_scores'] = seq_data['quality_scores'][start_pos:end_pos+1]
        
        if 'peak_positions' in seq_data and seq_data['peak_positions'] is not None:
            seq_data['peak_positions'] = seq_data['peak_positions'][start_pos:end_pos+1]
        
        # Trim trace data
        if 'plot_data' in seq_data:
            for base in seq_data['plot_data']:
                if len(seq_data['plot_data'][base]) > end_pos:
                    seq_data['plot_data'][base] = seq_data['plot_data'][base][start_pos:end_pos+1]
    
    def auto_correct_sequence(self, filename, confidence_threshold=0.8):
        """Automatically correct low-confidence base calls"""
        if filename not in self.app.sequences:
            return
        
        seq_data = self.app.sequences[filename]
        
        if 'quality_scores' not in seq_data or seq_data['quality_scores'] is None:
            messagebox.showwarning("Warning", "No quality scores available for auto-correction")
            return
        
        quality_scores = seq_data['quality_scores']
        sequence = seq_data['edited_seq']
        trace_data = seq_data.get('plot_data', {})
        
        corrections = 0
        
        for i, (base, quality) in enumerate(zip(sequence, quality_scores)):
            if quality < 20:  # Low quality threshold
                # Try to determine better base from trace data
                if trace_data:
                    best_base = self._get_best_base_from_traces(trace_data, i)
                    if best_base and best_base != base:
                        self.substitute_base(filename, i, best_base)
                        corrections += 1
        
        if corrections > 0:
            messagebox.showinfo("Auto-correction", f"Made {corrections} automatic corrections")
        else:
            messagebox.showinfo("Auto-correction", "No corrections needed")
    
    def _get_best_base_from_traces(self, trace_data, position):
        """Determine the best base call from trace data at a position"""
        if position >= min(len(trace) for trace in trace_data.values()):
            return None
        
        # Get signal intensities for each base
        intensities = {}
        for base, trace in trace_data.items():
            if position < len(trace):
                intensities[base] = trace[position]
        
        if not intensities:
            return None
        
        # Find base with highest intensity
        best_base = max(intensities, key=intensities.get)
        best_intensity = intensities[best_base]
        
        # Check if the signal is significantly higher than others
        other_intensities = [v for k, v in intensities.items() if k != best_base]
        if other_intensities:
            avg_other = np.mean(other_intensities)
            if best_intensity > avg_other * 2:  # At least 2x higher
                return best_base
        
        return None

