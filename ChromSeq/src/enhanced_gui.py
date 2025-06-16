"""
Enhanced GUI for Chromas Clone with Multi-Chromatogram Visualization
Provides ChromasPro-like functionality with interactive editing capabilities
"""

import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.widgets import SpanSelector, Cursor
import matplotlib.patches as patches
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import concurrent.futures
from collections import defaultdict
import threading
from tkinter import ttk
import pickle
from tkinter import filedialog, messagebox
# Import Bio.SeqIO for AB1 file parsing
from Bio import SeqIO

class EnhancedChromasClone:

    def __init__(self, root):
        self.root = root
        self.root.title("Enhanced Chromas Clone Pro")
        self.root.geometry("1600x1000")
        
        # Initialize data structures
        self.sequences = {}
        self.contigs = {}  # New: Store assembled contigs
        self.selected_sequences = []
        self.current_files = []
        self.zoom_level = 1.0
        self.pan_offset = 0
        self.alignment_result = None
        self.consensus = None
        self.project_file = None
        self.editing_mode = False
        self.selected_position = None
        
        # Enhanced visualization settings
        self.view_mode = "single"  # "single", "comparison", "overlay"
        self.show_quality = True
        self.show_peaks = True
        self.sync_navigation = True
        self.display_mode = "sequences"  # "sequences", "contigs", "both"
        
        # Base colors for display (enhanced with transparency)
        self.base_colors = {
            "A": "#00CC00",  # Green
            "C": "#0000FF",  # Blue
            "G": "#FF9900",  # Orange
            "T": "#FF0000",  # Red
            "N": "#808080"   # Gray
        }
        
        # Quality color mapping
        self.quality_colors = {
            "high": "#90EE90",    # Light green
            "medium": "#FFFF99",  # Light yellow
            "low": "#FFB6C1"      # Light pink
        }
        
        # Create executor for background tasks
        self.executor = concurrent.futures.ThreadPoolExecutor(max_workers=2)
        
        # Create status bar first, as other components might need to update it
        self.create_status_bar()

        # Initialize components after status bar is created
        from enhanced_utils import EnhancedUtils
        from enhanced_fileio import EnhancedFileIO
        from enhanced_sequence import EnhancedSequenceHandler
        from enhanced_alignment import EnhancedAlignment
        from enhanced_consensus import EnhancedConsensusGenerator
        
        self.utils = EnhancedUtils(self)
        self.file_io = EnhancedFileIO(self)
        self.sequence_handler = EnhancedSequenceHandler(self)
        self.alignment = EnhancedAlignment(self)
        self.consensus_gen = EnhancedConsensusGenerator(self)
        
        # Create GUI layout
        self.create_menu_bar()
        self.create_main_layout()
        self.create_toolbar()
        self.create_alignment_result_panel()

        
        # Bind events and shortcuts
        self.bind_shortcuts()
        self.setup_matplotlib_events()
        
        # Initialize view
        self.update_view_mode()


    def create_status_bar(self):
        """Create status bar at the bottom of the window"""
        self.status_bar = ttk.Label(self.root, text="Ready", relief=tk.SUNKEN, anchor=tk.W)
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)

        self.progress_var = tk.DoubleVar()
        self.progress_bar = ttk.Progressbar(
            self.status_bar,
            orient="horizontal",
            length=200,
            mode="determinate",
            variable=self.progress_var
        )
        self.progress_bar.pack(side=tk.RIGHT, padx=5)



    def load_project(self):
        """Load a saved Chromas Clone project (.chromasproj)"""
        file_path = filedialog.askopenfilename(
            title="Open Project",
            filetypes=[("Chromas Clone Project", "*.chromasproj"), ("All Files", "*.*")]
        )
        if not file_path:
            return

        try:
            with open(file_path, "rb") as f:
                project_data = pickle.load(f)

            # Restore main data structures
            self.app.sequences = project_data.get("sequences", {})
            self.app.contigs = project_data.get("contigs", {})
            self.app.project_file = file_path

            # Refresh GUI
            self.app.refresh_file_list()
            self.app.update_plot()
            messagebox.showinfo("Project Loaded", f"Successfully loaded: {file_path}")

        except Exception as e:
            messagebox.showerror("Load Error", f"Could not load project: {e}")

    def save_project(self):
        """Save current project to the existing file, or prompt Save As if none"""
        if not self.app.project_file:
            self.save_project_as()
            return

        try:
            project_data = {
                "sequences": self.app.sequences,
                "contigs": self.app.contigs,
            }
            with open(self.app.project_file, "wb") as f:
                pickle.dump(project_data, f)

            messagebox.showinfo("Project Saved", f"Project saved to: {self.app.project_file}")

        except Exception as e:
            messagebox.showerror("Save Error", f"Could not save project: {e}")

    def save_project_as(self):
        """Save current project to a new file"""
        file_path = filedialog.asksaveasfilename(
            title="Save Project As",
            defaultextension=".chromasproj",
            filetypes=[("Chromas Clone Project", "*.chromasproj"), ("All Files", "*.*")]
        )
        if not file_path:
            return

        self.app.project_file = file_path
        self.save_project()

    def export_fasta(self, selected_only=False, is_contig=False):
        """Export sequences or contigs to FASTA"""
        file_path = filedialog.asksaveasfilename(
            title="Export FASTA",
            defaultextension=".fasta",
            filetypes=[("FASTA Files", "*.fasta *.fa"), ("All Files", "*.*")]
        )
        if not file_path:
            return

        try:
            from Bio.SeqRecord import SeqRecord
            from Bio.Seq import Seq
            from Bio import SeqIO

            records = []
            if is_contig:
                for contig_name, contig_data in self.app.contigs.items():
                    seq = contig_data.get("consensus_sequence", "")
                    if seq:
                        records.append(SeqRecord(Seq(seq), id=contig_name, description="Contig Consensus"))
            else:
                if selected_only:
                    selected_items = self.app.file_tree.selection()
                    for item in selected_items:
                        name = self.app.file_tree.item(item, "text")
                        if name in self.app.sequences:
                            seq = self.app.sequences[name].get("edited_seq", "")
                            records.append(SeqRecord(Seq(seq), id=name, description="Edited Sequence"))
                else:
                    for name, data in self.app.sequences.items():
                        seq = data.get("edited_seq", "")
                        records.append(SeqRecord(Seq(seq), id=name, description="Edited Sequence"))

            SeqIO.write(records, file_path, "fasta")
            messagebox.showinfo("Export Complete", f"FASTA exported to: {file_path}")

        except Exception as e:
            messagebox.showerror("Export Error", f"Could not export FASTA: {e}")
    def create_menu_bar(self):
        """Create enhanced menu bar with all ChromasPro-like options"""
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        # File menu
        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Open AB1 Files...", command=self.file_io.add_files, accelerator="Ctrl+O")
        file_menu.add_command(label="Open Project...", command=self.load_project, accelerator="Ctrl+L")  # ✅ FIXED
        file_menu.add_separator()
        file_menu.add_command(label="Save Project", command=self.save_project, accelerator="Ctrl+S")    # ✅ FIXED
        file_menu.add_command(label="Save Project As...", command=self.save_project_as)                  # ✅ FIXED
        file_menu.add_separator()
        file_menu.add_command(label="Export FASTA...", command=self.export_fasta, accelerator="Ctrl+E") # ✅ FIXED
        file_menu.add_command(label="Export Image...", command=self.export_image)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)


        # Edit menu
        edit_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Edit", menu=edit_menu)
        edit_menu.add_command(label="Undo", command=self.sequence_handler.undo_edit, accelerator="Ctrl+Z")
        edit_menu.add_command(label="Redo", command=self.sequence_handler.redo_edit, accelerator="Ctrl+Y")
        edit_menu.add_separator()
        edit_menu.add_command(label="Toggle Edit Mode", command=self.toggle_editing_mode, accelerator="Ctrl+D")
        edit_menu.add_command(label="Insert Base...", command=self.sequence_handler.insert_base_dialog)
        edit_menu.add_command(label="Delete Base", command=self.sequence_handler.delete_selected_base)

        # View menu
        view_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="View", menu=view_menu)
        view_menu.add_command(label="Single View", command=lambda: self.set_view_mode("single"))
        view_menu.add_command(label="Comparison View", command=lambda: self.set_view_mode("comparison"))
        view_menu.add_command(label="Overlay View", command=lambda: self.set_view_mode("overlay"))
        view_menu.add_separator()
        view_menu.add_command(label="Show Sequences", command=lambda: self.set_display_mode("sequences"))
        view_menu.add_command(label="Show Contigs", command=lambda: self.set_display_mode("contigs"))
        view_menu.add_command(label="Show Both", command=lambda: self.set_display_mode("both"))
        view_menu.add_separator()
        view_menu.add_checkbutton(label="Show Quality", variable=tk.BooleanVar(value=self.show_quality),
                                  command=self.toggle_quality_display)
        view_menu.add_checkbutton(label="Show Peak Markers", variable=tk.BooleanVar(value=self.show_peaks),
                                  command=self.toggle_peak_display)
        view_menu.add_checkbutton(label="Sync Navigation", variable=tk.BooleanVar(value=self.sync_navigation),
                                  command=self.toggle_sync_navigation)

        # Analysis menu
        analysis_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Analysis", menu=analysis_menu)
        analysis_menu.add_command(label="Align Sequences", command=self.alignment.align_sequences, accelerator="Ctrl+A")
        analysis_menu.add_command(label="Generate Consensus", command=self.consensus_gen.generate_consensus)
        analysis_menu.add_separator()
        analysis_menu.add_command(label="Assemble All", command=self.assemble_all_sequences)
        analysis_menu.add_command(label="Assemble Selected", command=self.assemble_selected_sequences)
        analysis_menu.add_separator()
        analysis_menu.add_command(label="Quality Analysis", command=self.show_quality_analysis)
        analysis_menu.add_command(label="Trim Sequences", command=self.sequence_handler.trim_sequences)

        # Help menu
        help_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="User Guide", command=self.show_help)
        help_menu.add_command(label="About", command=self.show_about)

    
    def create_main_layout(self):
        """Create the main application layout with enhanced panels"""
        # Create main paned window
        main_paned = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        main_paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Left panel for file management and controls
        self.left_frame = ttk.Frame(main_paned, width=350)
        main_paned.add(self.left_frame, weight=0)
        
        # Right panel for visualization
        self.right_frame = ttk.Frame(main_paned)
        main_paned.add(self.right_frame, weight=1)
        
        # Create left panel components
        self.create_file_panel()
        self.create_control_panel()
        self.create_sequence_info_panel()
        
        # Create right panel components
        self.create_visualization_panel()
        self.create_consensus_panel()
    
    def create_file_panel(self):
        """Enhanced file management panel with role selection"""
        file_frame = ttk.LabelFrame(self.left_frame, text="Sequence Files & Contigs", padding="5")
        file_frame.pack(fill=tk.BOTH, padx=5, pady=5, expand=True)
        
        # Search and filter controls
        search_frame = ttk.Frame(file_frame)
        search_frame.pack(fill=tk.X, pady=(0,5))
        
        ttk.Label(search_frame, text="Search:").pack(side=tk.LEFT)
        self.search_var = tk.StringVar()
        self.search_var.trace("w", lambda *args: self.filter_files())
        search_entry = ttk.Entry(search_frame, textvariable=self.search_var)
        search_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(5,0))
        
        # Add selection info label
        self.selection_info = ttk.Label(search_frame, text="0 of 0 selected | Forward: 0 | Reverse: 0")
        self.selection_info.pack(side=tk.RIGHT, padx=(10,0))

        # Role assignment controls
        role_buttons_frame = ttk.Frame(file_frame)
        role_buttons_frame.pack(fill=tk.X, pady=(0,5))
        
        ttk.Button(role_buttons_frame, text="Forward",
                  command=lambda: self.assign_role_to_selected("forward")).pack(side=tk.LEFT, padx=(0,2))
        ttk.Button(role_buttons_frame, text="Reverse",
                  command=lambda: self.assign_role_to_selected("reverse")).pack(side=tk.LEFT, padx=(0,2))
        ttk.Button(role_buttons_frame, text="Clear Role",
                  command=lambda: self.assign_role_to_selected(None)).pack(side=tk.LEFT)
        
        # Assembly controls
        assembly_frame = ttk.Frame(file_frame)
        assembly_frame.pack(fill=tk.X, pady=(0,5))
        
        ttk.Button(assembly_frame, text="Assemble All",
                  command=self.assemble_all_sequences).pack(side=tk.LEFT, padx=(0,5))
        ttk.Button(assembly_frame, text="Assemble Selected",
                  command=self.assemble_selected_sequences).pack(side=tk.LEFT)
        
        # File list with enhanced display
        list_frame = ttk.Frame(file_frame)
        list_frame.pack(fill=tk.BOTH, expand=True)
        
        # Create Treeview for better file display
        columns = ("Name", "Type", "Role", "Length", "Quality")
        self.file_tree = ttk.Treeview(list_frame, columns=columns, show="tree headings", height=10)
        
        # Configure columns
        self.file_tree.heading("#0", text="File/Contig")
        self.file_tree.heading("Name", text="Name")
        self.file_tree.heading("Type", text="Type")
        self.file_tree.heading("Role", text="Role")
        self.file_tree.heading("Length", text="Length")
        self.file_tree.heading("Quality", text="Avg Quality")
        
        # Configure column widths
        self.file_tree.column("#0", width=150, minwidth=100)
        self.file_tree.column("Name", width=120, minwidth=80)
        self.file_tree.column("Type", width=60, minwidth=50)
        self.file_tree.column("Role", width=60, minwidth=50)
        self.file_tree.column("Length", width=80, minwidth=60)
        self.file_tree.column("Quality", width=80, minwidth=60)
        
        # Add scrollbars
        v_scrollbar = ttk.Scrollbar(list_frame, orient=tk.VERTICAL, command=self.file_tree.yview)
        h_scrollbar = ttk.Scrollbar(list_frame, orient=tk.HORIZONTAL, command=self.file_tree.xview)
        self.file_tree.configure(yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set)
        
        # Pack treeview and scrollbars
        self.file_tree.grid(row=0, column=0, sticky="nsew")
        v_scrollbar.grid(row=0, column=1, sticky="ns")
        h_scrollbar.grid(row=1, column=0, sticky="ew")
        
        list_frame.grid_rowconfigure(0, weight=1)
        list_frame.grid_columnconfigure(0, weight=1)
        
        # Bind events
        self.file_tree.bind("<<TreeviewSelect>>", self.on_file_selection)
        self.file_tree.bind("<Double-1>", self.on_file_double_click)
        self.file_tree.bind("<Button-3>", self.show_file_context_menu)  # Right-click
        
        # File management buttons
        button_frame = ttk.Frame(file_frame)
        button_frame.pack(fill=tk.X, pady=(5,0))
        
        ttk.Button(button_frame, text="Add Files", command=self.file_io.add_files).pack(side=tk.LEFT, padx=(0,5))
        ttk.Button(button_frame, text="Remove", command=self.remove_selected_files).pack(side=tk.LEFT, padx=(0,5))
        ttk.Button(button_frame, text="Refresh", command=self.refresh_file_list).pack(side=tk.LEFT)
    
    def assign_role_to_selected(self, role):
        """Assign role (forward/reverse/None) to selected sequences"""
        selected_items = self.file_tree.selection()
        if not selected_items:
            messagebox.showwarning("No Selection", "Please select sequences to assign roles.")
            return
        
        for item in selected_items:
            item_data = self.file_tree.item(item)
            filename = item_data["text"]
            
            # Only assign roles to sequences, not contigs
            if filename in self.sequences:
                self.sequences[filename]["role"] = role
                # Update display
                values = list(item_data["values"])
                values[2] = role if role else "None"  # Role column
                self.file_tree.item(item, values=values)
        
        self.utils.update_file_tree_colors()
        messagebox.showinfo("Role Assignment", f"Assigned role \'{role}\' to {len(selected_items)} sequences.")
    
    def assemble_all_sequences(self):
        """Assemble all sequences that have sufficient similarity"""
        if not self.sequences:
            messagebox.showwarning("No Sequences", "Please load sequences before assembly.")
            return
        
        # Start assembly in background
        self.utils.update_progress(0, 100, "Starting assembly of all sequences...")
        self.executor.submit(self._assemble_sequences_background, list(self.sequences.keys()))
    
    def assemble_selected_sequences(self):
        """Assemble selected sequences"""
        selected_items = self.file_tree.selection()
        if not selected_items:
            messagebox.showwarning("No Selection", "Please select sequences to assemble.")
            return
        
        # Get selected sequence names
        selected_names = []
        for item in selected_items:
            filename = self.file_tree.item(item, "text")
            if filename in self.sequences:  # Only sequences, not contigs
                selected_names.append(filename)
        
        if not selected_names:
            messagebox.showwarning("No Sequences", "Please select sequences (not contigs) to assemble.")
            return
        
        # Start assembly in background
        self.utils.update_progress(0, 100, f"Starting assembly of {len(selected_names)} selected sequences...")
        self.executor.submit(self._assemble_sequences_background, selected_names)
    
    def _assemble_sequences_background(self, sequence_names):
        """Background task for sequence assembly"""
        try:
            from enhanced_assembly import EnhancedAssembly
            assembler = EnhancedAssembly(self)
            
            # Perform assembly
            contigs = assembler.assemble_sequences(sequence_names)
            
            # Update GUI in main thread
            self.root.after(0, self._handle_assembly_results, contigs)
            
        except Exception as e:
            self.root.after(0, self._handle_assembly_error, str(e))
    
    def _handle_assembly_results(self, contigs):
        """Handle assembly results and update GUI"""
        if not self.root.winfo_exists():
            return

        if contigs:
            for contig_name, contig_data in contigs.items():
                self.contigs[contig_name] = contig_data
                # Move assembled sequences under the contig in the file tree
                for seq_name in contig_data["members"]:
                    if seq_name in self.sequences:
                        # Remove from root level
                        for item in self.file_tree.get_children():
                            if self.file_tree.item(item, "text") == seq_name:
                                self.file_tree.delete(item)
                                break
                self.refresh_file_list() # Re-add contigs and update sequence display
            messagebox.showinfo("Assembly Complete", f"Successfully assembled {len(contigs)} contig(s).")
        else:
            messagebox.showinfo("Assembly Complete", "No contigs were assembled from the selected sequences.")
        self.utils.reset_progress()
        self.utils.update_selection_info()

    def _handle_assembly_error(self, error_message):
        """Handle assembly errors and display to user"""
        if not self.root.winfo_exists():
            return
        messagebox.showerror("Assembly Error", f"An error occurred during assembly: {error_message}")
        self.utils.reset_progress()

    def create_control_panel(self):
        """Create control panel for visualization and editing"""
        control_frame = ttk.LabelFrame(self.left_frame, text="Controls", padding="5")
        control_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # View Mode selection
        view_mode_frame = ttk.Frame(control_frame)
        view_mode_frame.pack(fill=tk.X, pady=(0,5))
        ttk.Label(view_mode_frame, text="View Mode:").pack(side=tk.LEFT)
        
        self.view_mode_var = tk.StringVar(value=self.view_mode)
        ttk.Radiobutton(view_mode_frame, text="Single", variable=self.view_mode_var, value="single",
                        command=self.update_view_mode).pack(side=tk.LEFT, padx=(5,0))
        ttk.Radiobutton(view_mode_frame, text="Comparison", variable=self.view_mode_var, value="comparison",
                        command=self.update_view_mode).pack(side=tk.LEFT, padx=(5,0))
        ttk.Radiobutton(view_mode_frame, text="Overlay", variable=self.view_mode_var, value="overlay",
                        command=self.update_view_mode).pack(side=tk.LEFT, padx=(5,0))
        
        # Display Mode selection
        display_mode_frame = ttk.Frame(control_frame)
        display_mode_frame.pack(fill=tk.X, pady=(0,5))
        ttk.Label(display_mode_frame, text="Display:").pack(side=tk.LEFT)
        
        self.display_mode_var = tk.StringVar(value=self.display_mode)
        ttk.Radiobutton(display_mode_frame, text="Sequences", variable=self.display_mode_var, value="sequences",
                        command=self.update_display_mode).pack(side=tk.LEFT, padx=(5,0))
        ttk.Radiobutton(display_mode_frame, text="Contigs", variable=self.display_mode_var, value="contigs",
                        command=self.update_display_mode).pack(side=tk.LEFT, padx=(5,0))
        ttk.Radiobutton(display_mode_frame, text="Both", variable=self.display_mode_var, value="both",
                        command=self.update_display_mode).pack(side=tk.LEFT, padx=(5,0))
        
        # Zoom and Pan controls
        zoom_pan_frame = ttk.Frame(control_frame)
        zoom_pan_frame.pack(fill=tk.X, pady=(0,5))
        
        ttk.Label(zoom_pan_frame, text="Zoom:").pack(side=tk.LEFT)
        self.zoom_slider = ttk.Scale(zoom_pan_frame, from_=0.1, to_=5.0, orient=tk.HORIZONTAL,
                                     command=self.set_zoom_level, value=self.zoom_level)
        self.zoom_slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(5,0))
        
        ttk.Button(zoom_pan_frame, text="Reset Zoom", command=self.reset_zoom).pack(side=tk.LEFT, padx=(5,0))
        
        # Editing mode toggle
        edit_mode_frame = ttk.Frame(control_frame)
        edit_mode_frame.pack(fill=tk.X, pady=(0,5))
        
        self.edit_mode_var = tk.BooleanVar(value=self.editing_mode)
        ttk.Checkbutton(edit_mode_frame, text="Enable Editing Mode", variable=self.edit_mode_var,
                        command=self.toggle_editing_mode).pack(side=tk.LEFT)
        
        # Other controls (e.g., quality/peak display toggles)
        display_options_frame = ttk.Frame(control_frame)
        display_options_frame.pack(fill=tk.X, pady=(0,5))
        
        self.show_quality_var = tk.BooleanVar(value=self.show_quality)
        ttk.Checkbutton(display_options_frame, text="Show Quality", variable=self.show_quality_var,
                        command=self.toggle_quality_display).pack(side=tk.LEFT, padx=(0,10))
        
        self.show_peaks_var = tk.BooleanVar(value=self.show_peaks)
        ttk.Checkbutton(display_options_frame, text="Show Peaks", variable=self.show_peaks_var,
                        command=self.toggle_peak_display).pack(side=tk.LEFT)

    def create_sequence_info_panel(self):
        """Create panel for displaying sequence information"""
        info_paned = ttk.PanedWindow(self.left_frame, orient=tk.VERTICAL)
        info_paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # General Info
        general_info_frame = ttk.LabelFrame(info_paned, text="General Info", padding="5")
        info_paned.add(general_info_frame, weight=1)
        self.general_info = tk.Text(general_info_frame, wrap=tk.WORD, height=10, state=tk.DISABLED)
        self.general_info.pack(fill=tk.BOTH, expand=True)
        
        # Quality Info
        quality_info_frame = ttk.LabelFrame(info_paned, text="Quality Info", padding="5")
        info_paned.add(quality_info_frame, weight=1)
        self.quality_info = tk.Text(quality_info_frame, wrap=tk.WORD, height=10, state=tk.DISABLED)
        self.quality_info.pack(fill=tk.BOTH, expand=True)

    def _draw_alignment_chromatogram(self):
        """Draw chromatograms and aligned sequences directly in sequence_display."""
        if not self.current_alignment:
            return

        forward_aligned = self.current_alignment["forward_aligned"]
        reverse_aligned = self.current_alignment["reverse_aligned"]
        forward_data = self.current_alignment["forward_data"]
        reverse_data = self.current_alignment["reverse_data"]

        consensus_seq = ""
        match_line = ""

        for f, r in zip(forward_aligned, reverse_aligned):
            if f == r and f != "-":
                match_line += "|"
                consensus_seq += f
            elif f == "-" or r == "-":
                match_line += " "
                consensus_seq += "-"
            else:
                match_line += "."
                consensus_seq += "N"

        block_size = 50
        lines = []

        for i in range(0, len(forward_aligned), block_size):
            f_chunk = forward_aligned[i:i+block_size]
            m_chunk = match_line[i:i+block_size]
            r_chunk = reverse_aligned[i:i+block_size]
            c_chunk = consensus_seq[i:i+block_size]

            lines.append(f"Forward {i+1:<5} {f_chunk}")
            lines.append(f"         {' ' * 6} {m_chunk}")
            lines.append(f"Reverse {i+1:<5} {r_chunk}")
            lines.append(f"Consensus     {c_chunk}")
            lines.append("")

        # Update text display
        self.app.sequence_display.config(state=tk.NORMAL)
        self.app.sequence_display.delete("1.0", tk.END)
        self.app.sequence_display.insert(tk.END, "\n".join(lines))
        self.app.sequence_display.config(state=tk.DISABLED)

        # === Now plot chromatograms ===
        plot_data = {"forward": forward_data.get("plot_data", {}), "reverse": reverse_data.get("plot_data", {})}
        quality = {"forward": forward_data.get("quality_scores", []), "reverse": reverse_data.get("quality_scores", [])}
        base_colors = {"A": "green", "C": "blue", "G": "black", "T": "red"}

        self.app.axs[0].clear()
        self.app.axs[1].clear()
        self.app.axs[2].clear()

        for idx, role in enumerate(["forward", "reverse", "consensus"]):
            ax = self.app.axs[idx]
            ax.set_title(f"{role.capitalize()} trace")
            for base in "ACGT":
                if role == "consensus":
                    continue  # consensus chromatogram to be implemented if available
                signal = plot_data[role].get(base, [])
                ax.plot(signal, color=base_colors[base], label=base)
            ax.legend(loc="upper right")

        self.app.canvas.draw()
    def create_visualization_panel(self):
        """Create panel for chromatogram visualization"""
        viz_frame = ttk.LabelFrame(self.right_frame, text="Chromatogram View", padding="1")
        viz_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        self.fig, self.ax = plt.subplots(figsize=(6, 6))
        self.canvas = FigureCanvasTkAgg(self.fig, master=viz_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        self.toolbar = NavigationToolbar2Tk(self.canvas, viz_frame)
        self.toolbar.update()
        self.canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        # SpanSelector for region selection
        # rectprops argument is deprecated, remove it
        self.span = SpanSelector(self.ax, self.on_span_select, 'horizontal', useblit=True,
                                 button=1, interactive=True, drag_from_anywhere=True)
        
        # Cursor for base selection
        self.cursor = Cursor(self.ax, useblit=True, color='red', linewidth=1)
        
        # Text display for sequence and quality
        self.sequence_text_frame = ttk.Frame(viz_frame)
        self.sequence_text_frame.pack(fill=tk.X, pady=(5,0))
        
        self.sequence_display = tk.Text(self.sequence_text_frame, wrap=tk.NONE, height=5, font=("Courier New", 10))
        self.sequence_display.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        seq_h_scrollbar = ttk.Scrollbar(self.sequence_text_frame, orient=tk.HORIZONTAL, command=self.sequence_display.xview)
        seq_h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)
        self.sequence_display.config(xscrollcommand=seq_h_scrollbar.set)
        
        seq_v_scrollbar = ttk.Scrollbar(self.sequence_text_frame, orient=tk.VERTICAL, command=self.sequence_display.yview)
        seq_v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.sequence_display.config(yscrollcommand=seq_v_scrollbar.set)
    def create_consensus_panel(self):
        """Create panel for consensus sequence display"""
        consensus_frame = ttk.LabelFrame(self.right_frame, text="Consensus Sequence", padding="1")
        consensus_frame.pack(fill=tk.BOTH, expand=True, padx=15, pady=15)
        
        # Create the text widget and store it as an attribute
        self.consensus_display = tk.Text(
            consensus_frame, 
            wrap=tk.WORD, 
            height=5, 
            state=tk.DISABLED,
            font=('Courier New', 9)
        )
        
        # Add scrollbars
        scroll_y = ttk.Scrollbar(consensus_frame, orient=tk.VERTICAL, command=self.consensus_display.yview)
        scroll_x = ttk.Scrollbar(consensus_frame, orient=tk.HORIZONTAL, command=self.consensus_display.xview)
        self.consensus_display.configure(
            yscrollcommand=scroll_y.set,
            xscrollcommand=scroll_x.set
        )
        
        # Pack everything
        scroll_x.pack(side=tk.BOTTOM, fill=tk.X)
        scroll_y.pack(side=tk.RIGHT, fill=tk.Y)
        self.consensus_display.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
    def create_alignment_result_panel(self):
        """Create a panel for alignment result display"""
        alignment_frame = ttk.LabelFrame(self.right_frame, text="Alignment Result", padding="5")
        alignment_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Store as attribute so other modules can access it
        self.alignment_result = tk.Text(alignment_frame, wrap=tk.WORD, height=10, state=tk.DISABLED)
        scrollbar = ttk.Scrollbar(alignment_frame, orient=tk.VERTICAL, command=self.alignment_result.yview)
        self.alignment_result.configure(yscrollcommand=scrollbar.set)
        
        self.alignment_result.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

    def on_span_select(self, xmin, xmax):
        """Callback for region selection with SpanSelector."""
        start = int(round(min(xmin, xmax)))
        end = int(round(max(xmin, xmax)))
        max_len = self.get_max_sequence_length()
        start = max(0, min(start, max_len))
        end = max(0, min(end, max_len))

        self.selected_span = (start, end)

        # Clear previous span highlights if any
        if hasattr(self, 'span_patch'):
            self.span_patch.remove()

        self.span_patch = self.ax.axvspan(start, end, color='yellow', alpha=0.3)

        self.canvas.draw_idle()
        self.utils.update_status(f"Selected span: {start}–{end}")

    def create_toolbar(self):
        """Create a toolbar with quick access buttons"""
        toolbar_frame = ttk.Frame(self.root, relief=tk.RAISED)
        toolbar_frame.pack(side=tk.TOP, fill=tk.X)
        
        # Existing buttons
        ttk.Button(toolbar_frame, text="Open", command=self.file_io.add_files).pack(side=tk.LEFT, padx=2, pady=2)
        ttk.Button(toolbar_frame, text="Save", command=self.save_project).pack(side=tk.LEFT, padx=2, pady=2)
        ttk.Separator(toolbar_frame, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=2)
        
        # Alignment button
        ttk.Button(toolbar_frame, text="Align", command=self.alignment.align_sequences).pack(side=tk.LEFT, padx=2, pady=2)
        
        # Add the Consensus button here
        ttk.Button(toolbar_frame, 
                  text="Consensus", 
                  command=self.consensus_gen.generate_consensus).pack(side=tk.LEFT, padx=2, pady=2)
        
        ttk.Separator(toolbar_frame, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=2)
        
        # Remaining buttons
        ttk.Button(toolbar_frame, text="Zoom In", command=lambda: self.set_zoom_level(self.zoom_level * 1.2)).pack(side=tk.LEFT, padx=2, pady=2)
        ttk.Button(toolbar_frame, text="Zoom Out", command=lambda: self.set_zoom_level(self.zoom_level * 0.8)).pack(side=tk.LEFT, padx=2, pady=2)
        ttk.Button(toolbar_frame, text="Reset View", command=self.reset_zoom).pack(side=tk.LEFT, padx=2, pady=2)
        ttk.Separator(toolbar_frame, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=2)
        ttk.Button(toolbar_frame, text="Edit Mode", command=self.toggle_editing_mode).pack(side=tk.LEFT, padx=2, pady=2)

    def bind_shortcuts(self):
        """Bind keyboard shortcuts"""
        self.root.bind("<Control-o>", lambda event: self.file_io.add_files())
        self.root.bind("<Control-l>", lambda event: self.load_project())    # ✅ FIXED
        self.root.bind("<Control-s>", lambda event: self.save_project())    # ✅ FIXED
        self.root.bind("<Control-e>", lambda event: self.export_fasta())    # ✅ FIXED
        self.root.bind("<Control-z>", lambda event: self.sequence_handler.undo_edit())
        self.root.bind("<Control-y>", lambda event: self.sequence_handler.redo_edit())
        self.root.bind("<Control-d>", lambda event: self.toggle_editing_mode())
        self.root.bind("<Control-a>", lambda event: self.alignment.align_sequences())

    def setup_matplotlib_events(self):
        """Setup Matplotlib event handlers"""
        self.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.canvas.mpl_connect('button_press_event', self.on_mouse_click)
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)

    def on_scroll(self, event):
        """Handle scroll wheel for zooming"""
        if event.inaxes == self.ax:
            if event.button == 'up':
                self.set_zoom_level(self.zoom_level * 1.1)
            elif event.button == 'down':
                self.set_zoom_level(self.zoom_level * 0.9)

    def on_mouse_click(self, event):
        """Handle mouse clicks for base selection"""
        if event.inaxes == self.ax and self.editing_mode and event.button == 1:  # Left click
            x_data = int(round(event.xdata))
            self.selected_position = x_data
            self.draw_selection_marker()
            self.utils.update_status(f"Selected position: {self.selected_position}")
            
            # Update sequence display to highlight selected base
            self.update_sequence_display()

    def on_mouse_move(self, event):
        """Handle mouse movement for cursor display"""
        if event.inaxes == self.ax:
            x_data = int(round(event.xdata))
            self.utils.update_status(f"Current position: {x_data}")

    def draw_selection_marker(self):
        """Draw a vertical line at the selected position"""
        if hasattr(self, 'selection_line'):
            self.selection_line.remove()
        if self.selected_position is not None:
            self.selection_line = self.ax.axvline(self.selected_position, color='purple', linestyle='--', linewidth=1)
        self.canvas.draw_idle()

    def set_zoom_level(self, new_zoom):
        """Set the zoom level for the chromatogram"""
        self.zoom_level = max(0.1, min(new_zoom, 10.0))  # Limit zoom between 0.1x and 10x
        self.zoom_slider.set(self.zoom_level)
        self.update_plot()

    def reset_zoom(self):
        """Reset zoom and pan to default"""
        self.zoom_level = 1.0
        self.pan_offset = 0
        self.zoom_slider.set(self.zoom_level)
        self.ax.set_xlim(0, self.get_max_sequence_length())
        self.ax.set_ylim(bottom=0) # Reset y-axis to auto-scale
        self.update_plot()

    def toggle_editing_mode(self):
        """Toggle editing mode on/off"""
        self.editing_mode = self.edit_mode_var.get()
        if self.editing_mode:
            self.utils.update_status("Editing mode ENABLED. Click on chromatogram to select base.")
        else:
            self.utils.update_status("Editing mode DISABLED.")
            self.selected_position = None
            self.draw_selection_marker() # Remove marker
        self.update_plot()

    def toggle_quality_display(self):
        """Toggle display of quality scores"""
        self.show_quality = self.show_quality_var.get()
        self.update_plot()

    def toggle_peak_display(self):
        """Toggle display of peak markers"""
        self.show_peaks = self.show_peaks_var.get()
        self.update_plot()

    def toggle_sync_navigation(self):
        """Toggle synchronized navigation across multiple plots"""
        self.sync_navigation = not self.sync_navigation
        self.update_plot()

    def set_view_mode(self, mode):
        """Set the chromatogram view mode (single, comparison, overlay)"""
        self.view_mode = mode
        self.update_plot()

    def set_display_mode(self, mode):
        """Set the display mode for sequences/contigs"""
        self.display_mode = mode
        self.refresh_file_list()
        self.update_plot()

    def update_view_mode(self):
        """Update the plot based on the selected view mode"""
        self.view_mode = self.view_mode_var.get()
        self.update_plot()

    def update_display_mode(self):
        """Update the file tree and plot based on the selected display mode"""
        self.display_mode = self.display_mode_var.get()
        self.refresh_file_list()
        self.update_plot()

    def on_file_selection(self, event):
        """Handle file selection in the treeview"""
        selected_items = self.file_tree.selection()
        if selected_items:
            filename = self.file_tree.item(selected_items[0], "text")
            if filename in self.sequences:
                self.utils.update_sequence_info(filename)
                self.update_plot()
            elif filename in self.contigs: # If a contig is selected
                self.utils.update_sequence_info(filename) # Placeholder for contig info
                self.update_plot() # Plot contig consensus
        else:
            self.utils.clear_sequence_info()
            self.update_plot() # Clear plot
        self.utils.update_selection_info()

    def on_file_double_click(self, event):
        """Handle double-click on file in treeview (e.g., open in new window or set as primary)"""
        item = self.file_tree.selection()[0]
        filename = self.file_tree.item(item, "text")
        if filename in self.sequences:
            messagebox.showinfo("Double Click", f"Double clicked on sequence: {filename}")
            # Implement logic to set this as primary sequence for single view
        elif filename in self.contigs:
            messagebox.showinfo("Double Click", f"Double clicked on contig: {filename}")
            # Implement logic to display contig details

    def show_file_context_menu(self, event):
        """Show context menu for file tree items"""
        # Select item on right click
        item = self.file_tree.identify_row(event.y)
        if item:
            self.file_tree.selection_set(item)
            
            context_menu = tk.Menu(self.root, tearoff=0)
            item_data = self.file_tree.item(item)
            filename = item_data["text"]

            if filename in self.sequences: # Sequence specific options
                context_menu.add_command(label="Set as Forward", command=lambda: self.assign_role_to_selected("forward"))
                context_menu.add_command(label="Set as Reverse", command=lambda: self.assign_role_to_selected("reverse"))
                context_menu.add_command(label="Clear Role", command=lambda: self.assign_role_to_selected(None))
                context_menu.add_separator()
                context_menu.add_command(label="Remove Sequence", command=self.remove_selected_files)
                context_menu.add_command(label="Export FASTA", command=lambda: self.file_io.export_fasta(selected_only=True))
            elif filename in self.contigs: # Contig specific options
                context_menu.add_command(label="View Contig Details", command=lambda: messagebox.showinfo("Contig Details", f"Details for {filename}"))
                context_menu.add_command(label="Export Contig FASTA", command=lambda: self.file_io.export_fasta(selected_only=True, is_contig=True))
                context_menu.add_command(label="Remove Contig", command=self.remove_selected_files)

            try:
                context_menu.tk_popup(event.x_root, event.y_root)
            finally:
                context_menu.grab_release()

    def remove_selected_files(self):
        """Remove selected sequences or contigs from the application"""
        selected_items = self.file_tree.selection()
        if not selected_items:
            messagebox.showwarning("No Selection", "Please select files to remove.")
            return
        
        if not messagebox.askyesno("Confirm Removal", "Are you sure you want to remove the selected files?"):
            return
        
        for item in selected_items:
            filename = self.file_tree.item(item, "text")
            if filename in self.sequences:
                del self.sequences[filename]
            elif filename in self.contigs:
                del self.contigs[filename]
            self.file_tree.delete(item)
        
        self.utils.clear_sequence_info()
        self.update_plot()
        self.utils.update_selection_info()
        messagebox.showinfo("Removal Complete", f"Removed {len(selected_items)} file(s).")

    def refresh_file_list(self):
        """Refresh the file tree display based on current sequences and contigs"""
        if not self.root.winfo_exists():
            return

        # Clear existing items
        for item in self.file_tree.get_children():
            self.file_tree.delete(item)
        
        # Add sequences based on display mode
        if self.display_mode in ["sequences", "both"]:
            for filename, seq_data in self.sequences.items():
                # Only add if not part of an assembled contig that is also displayed
                is_part_of_displayed_contig = False
                if self.display_mode == "both":
                    for contig_data in self.contigs.values():
                        if filename in contig_data["members"]:
                            is_part_of_displayed_contig = True
                            break
                
                if not is_part_of_displayed_contig:
                    self._add_item_to_file_tree(filename, seq_data, "sequence")
        
        # Add contigs based on display mode
        if self.display_mode in ["contigs", "both"]:
            for contig_name, contig_data in self.contigs.items():
                contig_id = self._add_item_to_file_tree(contig_name, contig_data, "contig")
                # Add member sequences as children of the contig
                for member_seq_name in contig_data["members"]:
                    if member_seq_name in self.sequences:
                        member_seq_data = self.sequences[member_seq_name]
                        self._add_item_to_file_tree(member_seq_name, member_seq_data, "sequence", parent=contig_id)
        
        self.utils.update_file_tree_colors()
        self.utils.update_selection_info()

    def _add_item_to_file_tree(self, name, data, item_type, parent=""):
        """Helper to add an item (sequence or contig) to the file tree"""
        role = data.get("role", "None") if item_type == "sequence" else "N/A"
        length = len(data.get("edited_seq", "")) if item_type == "sequence" else len(data.get("consensus_sequence", ""))
        quality = f"{np.mean(data.get('quality_scores', [0])):.1f}" if item_type == "sequence" and data.get("quality_scores") is not None else "N/A"

        
        item_id = self.file_tree.insert(parent, tk.END, text=name, values=(name, item_type.capitalize(), role, length, quality))
        return item_id

    def update_plot(self):
        """Update the Matplotlib plot based on selected sequences/contigs and view mode"""
        self.ax.clear()
        selected_items = self.file_tree.selection()
        
        if not selected_items:
            self.ax.set_title("No Sequence Selected")
            self.canvas.draw_idle()
            self.update_sequence_display("")
            return
        
        primary_item_name = self.file_tree.item(selected_items[0], "text")
        primary_item_type = self.file_tree.item(selected_items[0], "values")[1].lower()

        if primary_item_type == "sequence":
            primary_data = self.sequences.get(primary_item_name)
            if not primary_data:
                self.ax.set_title("Error: Sequence data not found")
                self.canvas.draw_idle()
                self.update_sequence_display("")
                return
            self._plot_sequence(self.ax, primary_data, primary_item_name)
            self.ax.set_title(f"Chromatogram: {primary_item_name}")
            self.update_sequence_display(primary_data.get("edited_seq", ""))

        elif primary_item_type == "contig":
            contig_data = self.contigs.get(primary_item_name)
            if not contig_data:
                self.ax.set_title("Error: Contig data not found")
                self.canvas.draw_idle()
                self.update_sequence_display("")
                return
            
            # Plot consensus sequence of the contig
            if "consensus_sequence" in contig_data and "consensus_trace" in contig_data:
                self._plot_contig_consensus(self.ax, contig_data, primary_item_name)
                self.ax.set_title(f"Contig Consensus: {primary_item_name}")
                self.update_sequence_display(contig_data["consensus_sequence"])
            else:
                self.ax.set_title(f"Contig: {primary_item_name} (No Consensus Data)")
                self.update_sequence_display("")

        # Adjust x-axis limits based on zoom and pan
        max_len = self.get_max_sequence_length()
        current_xlim = self.ax.get_xlim()
        
        # Calculate new x-limits based on zoom_level and pan_offset
        # This needs to be more sophisticated for proper pan/zoom
        # For now, just apply zoom to the full length
        
        # Set x-axis limits
        self.ax.set_xlim(0, max_len * self.zoom_level)
        self.ax.set_ylim(bottom=0) # Ensure y-axis starts from 0
        
        self.canvas.draw_idle()

    def _plot_sequence(self, ax, seq_data, title):
        """Helper to plot a single sequence chromatogram"""
        if "plot_data" not in seq_data or not seq_data["plot_data"]:
            ax.text(0.5, 0.5, "No chromatogram data available",
                    horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
            return

        trace_data = seq_data["plot_data"]
        sequence = seq_data.get("edited_seq", "")
        quality_scores = seq_data.get("quality_scores")
        peak_locations = seq_data.get("peak_locations")

        x = np.arange(len(trace_data.get("A", [])))

        for base, color in self.base_colors.items():
            if base != "N": # Don't plot N trace
                ax.plot(x, trace_data.get(base, []), color=color, label=base, alpha=0.7)

        # Plot quality scores as background shading
        if self.show_quality and quality_scores is not None:
            for i, score in enumerate(quality_scores):
                color = self.utils.get_quality_color(score)
                ax.axvspan(i - 0.5, i + 0.5, facecolor=color, alpha=0.3, edgecolor='none')

        # Plot peak markers
        if self.show_peaks and peak_locations is not None:
            for i, loc in enumerate(peak_locations):
                if i < len(sequence):
                    base = sequence[i]
                    if base in trace_data and loc < len(trace_data[base]):
                        ax.plot(loc, trace_data[base][loc], 'o', color='black', markersize=4)
                        ax.text(loc, trace_data[base][loc] + 50, base, ha='center', va='bottom', fontsize=8)

        ax.set_xlabel("Base Position")
        ax.set_ylabel("Intensity")
        ax.legend()
        ax.grid(True)

    def _plot_contig_consensus(self, ax, contig_data, title):
        """Helper to plot a contig consensus chromatogram"""
        consensus_sequence = contig_data.get("consensus_sequence")
        consensus_trace = contig_data.get("consensus_trace")
        consensus_quality = contig_data.get("consensus_quality")
        
        if not consensus_trace or not consensus_sequence:
            ax.text(0.5, 0.5, "No consensus chromatogram data available",
                    horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
            return

        x = np.arange(len(consensus_trace.get("A", [])))

        for base, color in self.base_colors.items():
            if base != "N":
                ax.plot(x, consensus_trace.get(base, []), color=color, label=base, alpha=0.7)

        if self.show_quality and consensus_quality is not None:
            for i, score in enumerate(consensus_quality):
                color = self.utils.get_quality_color(score)
                ax.axvspan(i - 0.5, i + 0.5, facecolor=color, alpha=0.3, edgecolor='none')

        # Add base calls to plot
        for i, base in enumerate(consensus_sequence):
            if base in consensus_trace and i < len(consensus_trace[base]):
                ax.text(i, consensus_trace[base][i] + 50, base, ha='center', va='bottom', fontsize=8)

        ax.set_xlabel("Base Position")
        ax.set_ylabel("Intensity")
        ax.legend()
        ax.grid(True)

    def update_sequence_display(self, sequence_str=""):
        """Update the sequence display text widget"""
        if not self.root.winfo_exists():
            return

        self.sequence_display.config(state=tk.NORMAL)
        self.sequence_display.delete(1.0, tk.END)
        
        if not sequence_str:
            self.sequence_display.config(state=tk.DISABLED)
            return

        formatted_seq = self.utils.format_sequence_for_display(sequence_str)
        self.sequence_display.insert(tk.END, formatted_seq)
        
        # Highlight selected base if in editing mode
        if self.editing_mode and self.selected_position is not None:
            line_length = 50 # Must match format_sequence_for_display
            line_num = self.selected_position // line_length
            char_in_line = self.selected_position % line_length
            
            # Calculate Tkinter text index
            # Each line has 6 chars for line number + ": " + line_seq + "\n"
            # So, line_num + 1 for 1-based line index
            # 8 for the prefix "   123: "
            start_index = f"{line_num + 1}.{8 + char_in_line}"
            end_index = f"{line_num + 1}.{8 + char_in_line + 1}"
            
            self.sequence_display.tag_add("highlight", start_index, end_index)
            self.sequence_display.tag_config("highlight", background="yellow")
            
            # Scroll to selected position
            self.sequence_display.see(start_index)

        self.sequence_display.config(state=tk.DISABLED)

    def get_max_sequence_length(self):
        """Get the maximum length among currently displayed sequences/contigs"""
        max_len = 0
        if self.display_mode in ["sequences", "both"]:
            for seq_data in self.sequences.values():
                max_len = max(max_len, len(seq_data.get("edited_seq", "")))
        if self.display_mode in ["contigs", "both"]:
            for contig_data in self.contigs.values():
                max_len = max(max_len, len(contig_data.get("consensus_sequence", "")))
        return max_len

    def export_image(self):
        """Export the current chromatogram view as an image file"""
        file_path = filedialog.asksaveasfilename(defaultextension=".png",
                                               filetypes=[("PNG files", "*.png"),
                                                          ("JPEG files", "*.jpg"),
                                                          ("All files", "*.*")],
                                               title="Export Chromatogram Image")
        if file_path:
            try:
                self.fig.savefig(file_path, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Export Complete", f"Chromatogram exported to {file_path}")
            except Exception as e:
                messagebox.showerror("Export Error", f"Failed to export image: {e}")

    def show_quality_analysis(self):
        """Display a detailed quality analysis report"""
        messagebox.showinfo("Quality Analysis", "Detailed quality analysis report will be displayed here.")

    def show_help(self):
        """Display help documentation"""
        messagebox.showinfo("Help", "User guide and help information.")

    def show_about(self):
        """Display about information"""
        messagebox.showinfo("About", "Enhanced Chromas Clone Pro\nVersion 1.0\nDeveloped by Manus")


if __name__ == "__main__":
    root = tk.Tk()
    app = EnhancedChromasClone(root)
    root.mainloop()


