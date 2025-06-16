"""
Enhanced utilities for the improved Chromas Clone application
Provides advanced file management and UI update functionality
"""

import tkinter as tk
from tkinter import ttk
import numpy as np

class EnhancedUtils:
    def __init__(self, app):
        self.app = app

    def get_quality_color(self, score):
        """
        Return a hex color based on a Phred quality score.
        Higher = greener, lower = redder.
        """
        if score >= 30:
            return "#00CC00"  # green
        elif score >= 20:
            return "#CCCC00"  # yellow
        else:
            return "#CC0000"  # red

    def clear_sequence_info(self):
        """Clear both general and quality info text boxes."""
        if not self.app.root.winfo_exists():
            return

        self.app.general_info.config(state=tk.NORMAL)
        self.app.general_info.delete(1.0, tk.END)
        self.app.general_info.config(state=tk.DISABLED)

        self.app.quality_info.config(state=tk.NORMAL)
        self.app.quality_info.delete(1.0, tk.END)
        self.app.quality_info.config(state=tk.DISABLED)

    def format_sequence_for_display(self, sequence_str, line_length=80):
        """
        Split a long sequence string into lines of fixed length for better display.
        """
        return "\n".join(sequence_str[i:i+line_length] for i in range(0, len(sequence_str), line_length))


    def update_selection_info(self):
        """Update selection information display with enhanced details"""
        selected_items = self.app.file_tree.selection()
        total_files = len(self.app.sequences)
        
        # Count sequences by role
        forward_count = sum(1 for s in self.app.sequences.values() if s.get("role") == "forward")
        reverse_count = sum(1 for s in self.app.sequences.values() if s.get("role") == "reverse")
        
        # Calculate total length and average quality
        total_length = sum(len(s.get("edited_seq", "")) for s in self.app.sequences.values())
        avg_quality = self.calculate_average_quality()
        
        # Update display
        info_text = (f"{len(selected_items)} of {total_files} selected | " 
                     f"Forward: {forward_count} | Reverse: {reverse_count}")
        if total_length > 0:
            info_text += f" | Total: {total_length} bp"
        if avg_quality is not None:
            info_text += f" | Avg Quality: {avg_quality:.1f}"
        
        if self.app.root.winfo_exists():
            self.app.selection_info.config(text=info_text)
        
        # Update file tree display with role coloring
        self.update_file_tree_colors()
    
    def calculate_average_quality(self):
        """Calculate average quality score across all sequences"""
        if not self.app.sequences:
            return None
        
        total_quality = 0
        total_bases = 0
        
        for seq_data in self.app.sequences.values():
            if "quality_scores" in seq_data and seq_data["quality_scores"] is not None:
                quality_scores = seq_data["quality_scores"]
                total_quality += np.sum(quality_scores)
                total_bases += len(quality_scores)
        
        return total_quality / total_bases if total_bases > 0 else None
    
    def update_file_tree_colors(self):
        """Update file tree with role-based coloring"""
        if not self.app.root.winfo_exists():
            return

        for item in self.app.file_tree.get_children():
            filename = self.app.file_tree.item(item, "text")
            if filename in self.app.sequences:
                role = self.app.sequences[filename].get("role")
                
                # Set tags based on role
                if role == "forward":
                    self.app.file_tree.item(item, tags=("forward",))
                elif role == "reverse":
                    self.app.file_tree.item(item, tags=("reverse",))
                else:
                    self.app.file_tree.item(item, tags=("normal",))
        
        # Configure tag colors
        self.app.file_tree.tag_configure("forward", background="#E8F5E8")
        self.app.file_tree.tag_configure("reverse", background="#FFE8E8")
        self.app.file_tree.tag_configure("normal", background="white")
    
    def update_sequence_info(self, filename=None):
        """Update sequence information panels with detailed data"""
        if not self.app.root.winfo_exists():
            return

        if filename is None:
            selected_items = self.app.file_tree.selection()
            if not selected_items:
                self.clear_sequence_info()
                return
            filename = self.app.file_tree.item(selected_items[0], "text")
        
        if filename not in self.app.sequences:
            return
        
        seq_data = self.app.sequences[filename]
        
        # Update general info
        self.update_general_info(filename, seq_data)
        
        # Update quality info
        self.update_quality_info(filename, seq_data)
    def update_general_info(self, filename, seq_data):
        """Update general sequence information"""
        if not self.app.root.winfo_exists():
            return

        self.app.general_info.config(state=tk.NORMAL)
        self.app.general_info.delete(1.0, tk.END)

        info_text = f"File: {filename}\n"
        info_text += f"Role: {seq_data.get('role', 'None')}\n"

        if "seq_record" in seq_data:
            seq_record = seq_data["seq_record"]
            info_text += f"Original Length: {len(seq_record.seq)} bp\n"

            # Add ABIF metadata if available
            if hasattr(seq_record, "annotations") and "abif_raw" in seq_record.annotations:
                abif_data = seq_record.annotations["abif_raw"]

                # Sample name
                if "SMPL1" in abif_data:
                    info_text += f"Sample: {abif_data['SMPL1']}\n"

                # Run date/time
                if "RUND1" in abif_data:
                    info_text += f"Run Date: {abif_data['RUND1']}\n"
                if "RUNT1" in abif_data:
                    info_text += f"Run Time: {abif_data['RUNT1']}\n"

                # Machine info
                if "MCHN1" in abif_data:
                    info_text += f"Machine: {abif_data['MCHN1']}\n"

                # Mobility file
                if "MOBF1" in abif_data:
                    info_text += f"Mobility File: {abif_data['MOBF1']}\n"

        # Base composition
        sequence = seq_data.get("edited_seq", "")
        if sequence:
            base_counts = {base: sequence.count(base) for base in "ATCG"}
            total_bases = sum(base_counts.values())

            info_text += "\nBase Composition:\n"
            for base in "ATCG":
                count = base_counts[base]
                percentage = (count / total_bases * 100) if total_bases > 0 else 0
                info_text += f"  {base}: {count} ({percentage:.1f}%)\n"

        self.app.general_info.insert(1.0, info_text)
        self.app.general_info.config(state=tk.DISABLED)

    
    def update_quality_info(self, filename, seq_data):
        """Update quality information"""
        if not self.app.root.winfo_exists():
            return

        self.app.quality_info.config(state=tk.NORMAL)
        self.app.quality_info.delete(1.0, tk.END)
        
        quality_text = f"Quality Analysis for {filename}\n\n"
        
        # Extract quality scores if available
        quality_scores = seq_data.get("quality_scores")
        if quality_scores is not None:
            avg_quality = np.mean(quality_scores)
            min_quality = np.min(quality_scores)
            max_quality = np.max(quality_scores)
            
            quality_text += f"Average Quality: {avg_quality:.2f}\n"
            quality_text += f"Min Quality: {min_quality:.2f}\n"
            quality_text += f"Max Quality: {max_quality:.2f}\n\n"
            
            # Quality distribution
            high_quality = np.sum(quality_scores >= 30)
            medium_quality = np.sum((quality_scores >= 20) & (quality_scores < 30))
            low_quality = np.sum(quality_scores < 20)
            total = len(quality_scores)
            
            quality_text += "Quality Distribution:\n"
            quality_text += f"  High (≥30): {high_quality} ({high_quality/total*100:.1f}%)\n"
            quality_text += f"  Medium (20-29): {medium_quality} ({medium_quality/total*100:.1f}%)\n"
            quality_text += f"  Low (<20): {low_quality} ({low_quality/total*100:.1f}%)\n\n"
            
            # Identify problematic regions
            low_quality_regions = self.find_low_quality_regions(quality_scores)
            if low_quality_regions:
                quality_text += "Low Quality Regions:\n"
                for start, end in low_quality_regions:
                    quality_text += f"  {start}-{end} (length: {end-start+1})\n"
        else:
            quality_text += "No quality scores available\n"
        
        # Peak analysis if trace data available
        if "plot_data" in seq_data:
            peak_analysis = self.analyze_peaks(seq_data["plot_data"])
            quality_text += "\nPeak Analysis:\n"
            quality_text += peak_analysis
        
        self.app.quality_info.insert(1.0, quality_text)
        self.app.quality_info.config(state=tk.DISABLED)
    
    def find_low_quality_regions(self, quality_scores, threshold=20, min_length=5):
        """Find regions with consistently low quality scores"""
        low_quality_regions = []
        in_region = False
        region_start = 0
        
        for i, score in enumerate(quality_scores):
            if score < threshold:
                if not in_region:
                    in_region = True
                    region_start = i
            else:
                if in_region:
                    region_length = i - region_start
                    if region_length >= min_length:
                        low_quality_regions.append((region_start, i-1))
                    in_region = False
        
        # Handle case where sequence ends in low quality region
        if in_region:
            region_length = len(quality_scores) - region_start
            if region_length >= min_length:
                low_quality_regions.append((region_start, len(quality_scores)-1))
        
        return low_quality_regions
    
    def analyze_peaks(self, trace_data):
        """Analyze peak characteristics in trace data"""
        analysis = ""
        
        for base, values in trace_data.items():
            if len(values) == 0:
                continue
                
            max_intensity = np.max(values)
            avg_intensity = np.mean(values)
            peak_count = self.count_peaks(values)
            
            analysis += f"  {base}: Max={max_intensity:.0f}, Avg={avg_intensity:.0f}, Peaks={peak_count}\n"
        
        return analysis
    
    def count_peaks(self, values, min_height=None, min_distance=1):
        """Count peaks in signal data"""
        if min_height is None:
            min_height = np.mean(values) + np.std(values)
        
        peaks = 0
        in_peak = False
        last_peak_pos = -min_distance - 1
        
        for i, value in enumerate(values):
            if value >= min_height and i - last_peak_pos > min_distance:
                if not in_peak:
                    peaks += 1
                    in_peak = True
                    last_peak_pos = i
            elif value < min_height * 0.8:  # Peak end threshold
                in_peak = False
        
        return peaks
    
    def update_progress(self, value, maximum=100, text=""):
        """Update progress bar and status"""
        if self.app.root.winfo_exists():
            # Schedule the GUI update to run on the main thread
            self.app.root.after(0, self._update_progress_gui, value, maximum, text)
    
    def _update_progress_gui(self, value, maximum, text):
        """Internal method to update progress bar and status on the main thread."""
        # Check again if root exists, as it might have been destroyed between after() call and execution
        if self.app.root.winfo_exists() and hasattr(self.app, "status_bar") and self.app.status_bar.winfo_exists():
            self.app.progress_var.set(value)
            if text:
                self.app.status_bar.config(text=text)

    def reset_progress(self):
        """Reset progress bar"""
        if self.app.root.winfo_exists():
            # Schedule the GUI update to run on the main thread
            self.app.root.after(0, self._reset_progress_gui)

    def _reset_progress_gui(self):
        """Internal method to reset progress bar on the main thread."""
        if self.app.root.winfo_exists() and hasattr(self.app, "status_bar") and self.app.status_bar.winfo_exists():
            self.app.progress_var.set(0)
            self.app.status_bar.config(text="Ready")

    def update_status(self, text):
        """Update the status bar text"""
        if self.app.root.winfo_exists():
            # Schedule the GUI update to run on the main thread
            self.app.root.after(0, lambda: self._update_status_gui(text))

    def _update_status_gui(self, text):
        """Internal method to update status bar on the main thread."""
        if self.app.root.winfo_exists() and hasattr(self.app, "status_bar") and self.app.status_bar.winfo_exists():
            self.app.status_bar.config(text=text)




