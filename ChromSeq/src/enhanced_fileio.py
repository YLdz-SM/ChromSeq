"""
Enhanced file I/O operations for the improved Chromas Clone application.
Supports advanced AB1 parsing, project management, and export functionality.
"""

import os
import pickle
import struct
import numpy as np
from tkinter import filedialog, messagebox
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from abifpy import Trace  # ✅ Correct import

class EnhancedFileIO:
    def __init__(self, app):
        self.app = app
        self.compression_enabled = True
        self.supported_formats = {
            'ab1': 'ABI Sequencer Files',
            'scf': 'SCF Files',
            'fasta': 'FASTA Files',
            'fastq': 'FASTQ Files'
        }

    def add_files(self):
        """Open file dialog and load files in background"""
        filetypes = [
            ("ABI Files", "*.ab1"),
            ("All Files", "*")
        ]

        try:
            self.app.root.update_idletasks()
            filepaths = filedialog.askopenfilenames(
                title="Select Sequence Files",
                filetypes=filetypes
            )
            self.app.root.update_idletasks()
        except Exception as e:
            messagebox.showerror("File Selection Error", f"Error: {e}")
            return

        if not filepaths:
            return

        self.app.root.after(0, self.app.utils.update_progress, 0, len(filepaths), "Loading files...")
        self.app.executor.submit(self._load_files_background, filepaths)

    def _load_files_background(self, filepaths):
        """Load multiple files in a thread"""
        loaded = 0
        failed = []

        for i, fp in enumerate(filepaths):
            if not self.app.root.winfo_exists():
                break
            try:
                fname, fdata = self._load_single_file(fp)
                if fname and fdata:
                    self.app.root.after(0, self._add_loaded_file, fname, fdata)
                    loaded += 1
                else:
                    failed.append(os.path.basename(fp))
            except Exception as e:
                failed.append(f"{os.path.basename(fp)}: {e}")

            if self.app.root.winfo_exists():
                progress = ((i + 1) / len(filepaths)) * 100
                self.app.root.after(0, self.app.utils.update_progress, progress, 100, f"Loading... {i+1}/{len(filepaths)}")

        if self.app.root.winfo_exists():
            self.app.root.after(0, self._show_load_results, loaded, failed)
            self.app.root.after(0, self.app.utils.reset_progress)

    def _add_loaded_file(self, filename, data):
        """Add one loaded file to app storage"""
        self.app.sequences[filename] = data
        self.app.refresh_file_list()
        self.app.update_plot()

    def _show_load_results(self, loaded, failed):
        """Show result message"""
        msg = f"Loaded: {loaded} file(s).\n"
        if failed:
            msg += f"Failed: {len(failed)}\n" + "\n".join(failed)
        messagebox.showinfo("Load Results", msg)

    def _load_single_file(self, filepath):
        """Detect and load one file"""
        ext = os.path.splitext(filepath)[1].lower()
        if ext == '.ab1':
            return self._load_ab1_file(filepath)
        elif ext in ['.fasta', '.fas', '.fa']:
            return self._load_fasta_file(filepath)
        elif ext in ['.fastq', '.fq']:
            return self._load_fastq_file(filepath)
        else:
            return None, None

    def _load_ab1_file(self, filepath):
        """Load AB1 using abifpy.Trace with fallback"""
        filename = os.path.basename(filepath)
        try:
            abif = Trace(filepath)
            seq = abif.seq
            length = len(seq)

            trace_data = {}
            for base, tag in zip("GATC", [9, 10, 11, 12]):
                raw = abif.data.get(('DATA', tag), b'')
                trace_data[base] = np.frombuffer(raw, dtype='>h')

            ploc = abif.data.get(('PLOC', 2), b'')
            peak_positions = np.frombuffer(ploc, dtype='>h')[:length]

            qual = np.array(abif.qual_val, dtype=np.uint8) if hasattr(abif, 'qual_val') else np.zeros(length, dtype=np.uint8)
            metadata = {f"{k[0]}{k[1]}": abif.data[k] for k in abif.data.keys()}

            return filename, {
                'seq_record': SeqRecord(Seq(seq), id=filename),
                'edited_seq': seq,
                'original_seq': seq,
                'plot_data': trace_data,
                'quality_scores': qual,
                'peak_positions': peak_positions,
                'metadata': metadata,
                'filepath': filepath,
                'role': None,
                'actual_length': length,
                'file_type': 'ab1',
            }
        except Exception as e:
            print(f"[WARN] abifpy failed: {e}")
            try:
                rec = SeqIO.read(filepath, "abi")
                seq = str(rec.seq)
                raw = rec.annotations.get('abif_raw', {})
                trace_data = {}
                for base, tag in zip("GATC", [9, 10, 11, 12]):
                    arr = raw.get(f'DATA{tag}', b'')
                    trace_data[base] = np.frombuffer(arr, dtype='>h')
                ploc = raw.get('PLOC2', b'')
                peak_positions = np.frombuffer(ploc, dtype='>h')[:len(seq)]
                qual = np.array(rec.letter_annotations.get('phred_quality', []), dtype=np.uint8)
                return filename, {
                    'seq_record': rec,
                    'edited_seq': seq,
                    'original_seq': seq,
                    'plot_data': trace_data,
                    'quality_scores': qual,
                    'peak_positions': peak_positions,
                    'metadata': raw,
                    'filepath': filepath,
                    'role': None,
                    'actual_length': len(seq),
                    'file_type': 'ab1',
                }
            except Exception as e2:
                print(f"[ERROR] SeqIO fallback failed: {e2}")
                return None, None

    def _load_fasta_file(self, filepath):
        """Simple FASTA loader"""
        recs = list(SeqIO.parse(filepath, "fasta"))
        if recs:
            rec = recs[0]
            seq = str(rec.seq)
            return os.path.basename(filepath), {
                'seq_record': rec,
                'edited_seq': seq,
                'original_seq': seq,
                'plot_data': {},
                'quality_scores': np.zeros(len(seq), dtype=np.uint8),
                'peak_positions': np.arange(len(seq)),
                'metadata': {},
                'filepath': filepath,
                'role': None,
                'actual_length': len(seq),
                'file_type': 'fasta',
            }
        return None, None

    def _load_fastq_file(self, filepath):
        """Simple FASTQ loader"""
        recs = list(SeqIO.parse(filepath, "fastq"))
        if recs:
            rec = recs[0]
            seq = str(rec.seq)
            qual = np.array(rec.letter_annotations.get('phred_quality', []), dtype=np.uint8)
            return os.path.basename(filepath), {
                'seq_record': rec,
                'edited_seq': seq,
                'original_seq': seq,
                'plot_data': {},
                'quality_scores': qual,
                'peak_positions': np.arange(len(seq)),
                'metadata': {},
                'filepath': filepath,
                'role': None,
                'actual_length': len(seq),
                'file_type': 'fastq',
            }
        return None, None

    def save_project(self):
        """Save project"""
        if not self.app.project_file:
            self.save_project_as()
            return
        try:
            with open(self.app.project_file, "wb") as f:
                pickle.dump({
                    'sequences': self.app.sequences,
                    'contigs': self.app.contigs,
                }, f)
            messagebox.showinfo("Project Saved", f"Saved to: {self.app.project_file}")
        except Exception as e:
            messagebox.showerror("Save Error", f"Could not save project: {e}")

    def save_project_as(self):
        """Save project with new filename"""
        fp = filedialog.asksaveasfilename(
            title="Save Project As",
            defaultextension=".chromasproj",
            filetypes=[("Chromas Clone Project", "*.chromasproj"), ("All Files", "*.*")]
        )
        if not fp:
            return
        self.app.project_file = fp
        self.save_project()

    def export_fasta(self, selected_only=False, is_contig=False):
        """Export sequences to FASTA"""
        fp = filedialog.asksaveasfilename(
            title="Export FASTA",
            defaultextension=".fasta",
            filetypes=[("FASTA Files", "*.fasta *.fa"), ("All Files", "*.*")]
        )
        if not fp:
            return

        try:
            records = []
            if is_contig:
                for cname, cdata in self.app.contigs.items():
                    seq = cdata.get('consensus_sequence', '')
                    if seq:
                        records.append(SeqRecord(Seq(seq), id=cname))
            else:
                if selected_only:
                    sel = self.app.file_tree.selection()
                    for item in sel:
                        name = self.app.file_tree.item(item, "text")
                        if name in self.app.sequences:
                            seq = self.app.sequences[name]['edited_seq']
                            records.append(SeqRecord(Seq(seq), id=name))
                else:
                    for name, data in self.app.sequences.items():
                        seq = data['edited_seq']
                        records.append(SeqRecord(Seq(seq), id=name))
            SeqIO.write(records, fp, "fasta")
            messagebox.showinfo("Export Complete", f"FASTA saved to: {fp}")
        except Exception as e:
            messagebox.showerror("Export Error", f"Could not export: {e}")
