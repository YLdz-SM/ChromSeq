#!/usr/bin/env python3
"""
Enhanced Chromas Clone Pro - Complete Application
Main entry point for the enhanced sequence analysis application
"""

import sys
import os
import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib
matplotlib.use('TkAgg')

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

try:
    from enhanced_gui import EnhancedChromasClone
except ImportError as e:
    print(f"Error importing enhanced GUI: {e}")
    print("Please ensure all required packages are installed:")
    print("pip install -r requirements.txt")
    sys.exit(1)

def main():
    """Main application entry point"""
    try:
        # Create and configure root window
        root = tk.Tk()
        root.title("Enhanced Chromas Clone Pro")
        root.geometry("1200x800")
        
        # Set application icon (if available)
        try:
            # You can add an icon file here
            # root.iconbitmap('icon.ico')
            pass
        except:
            pass
        
        # Create application instance
        app = EnhancedChromasClone(root)
        
        # Configure window closing
        def on_closing():
            if app.sequences:
                if messagebox.askyesno("Quit", "Do you want to save your project before closing?"):
                    app.file_io.save_project()
            root.destroy()
        
        root.protocol("WM_DELETE_WINDOW", on_closing)
        
        # Start the application
        print("Enhanced Chromas Clone Pro starting...")
        print("Loading complete. Ready for sequence analysis.")
        
        root.mainloop()
        
    except Exception as e:
        error_msg = f"Failed to start Enhanced Chromas Clone Pro:\\n{str(e)}"
        print(error_msg)
        
        # Try to show error in GUI if possible
        try:
            root = tk.Tk()
            root.withdraw()
            messagebox.showerror("Startup Error", error_msg)
            root.destroy()
        except:
            pass
        
        sys.exit(1)

if __name__ == "__main__":
    main()

