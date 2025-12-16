#!/usr/bin/env python3
"""
Hybrid Metabolism Prediction GUI

Combines two approaches:
1. Simple pattern matching (your current metabolic_liability_tagger)
2. SMARTCyp-inspired quantitative scoring

Features:
- Side-by-side comparison
- Method toggle
- Quantitative ranking
- Dark theme UI
"""

import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext
import threading
import sys
import io
from pathlib import Path


class HybridMetabolismGUI:
    """
    Hybrid GUI for comparing metabolism prediction methods
    """

    # Dark theme colors
    COLORS = {
        'bg_dark': '#1a1a1a',
        'bg_medium': '#2d2d2d',
        'bg_light': '#3d3d3d',
        'input_bg': '#ffffff',
        'input_fg': '#000000',
        'text_primary': '#ffffff',
        'text_secondary': '#b0b0b0',
        'text_on_button': '#000000',  # Black text for buttons
        'btn_primary': '#4CAF50',     # Green primary button
        'btn_primary_hover': '#45a049',
        'btn_success': '#2196F3',     # Blue for success
        'btn_info': '#2196F3',
        'btn_warning': '#ff9800',     # Orange for warning
        'accent_green': '#4CAF50',
        'accent_blue': '#2196F3',
        'accent_orange': '#ff9800',
        'accent_purple': '#9C27B0',
    }

    def __init__(self, root):
        self.root = root
        self.root.title('Hybrid Metabolism Predictor - Pattern vs Quantitative')
        self.root.geometry('700x750')
        self.root.resizable(False, False)
        self.root.configure(bg=self.COLORS['bg_dark'])

        # Default test SMILES
        self.default_smiles = "CC(C)[C@@H](C(=O)N1Cc2ccccc2[C@H]1C(=O)Nc1cn(C)c(=O)n(C)c1=O)NC(=O)[C@H](CC1CCCCC1)NC(=O)C"

        # Prediction method selection
        self.method_var = tk.StringVar(value="quantitative")

        self.setup_ui()
        self.center_window()

    def setup_ui(self):
        """Create UI elements"""

        # Main container
        main_frame = tk.Frame(
            self.root,
            bg=self.COLORS['bg_dark'],
            padx=25,
            pady=20
        )
        main_frame.pack(fill=tk.BOTH, expand=True)

        # ========== HEADER ==========
        header_frame = tk.Frame(main_frame, bg=self.COLORS['bg_dark'])
        header_frame.pack(fill=tk.X, pady=(0, 15))

        title_label = tk.Label(
            header_frame,
            text="HYBRID METABOLISM PREDICTOR",
            font=('Helvetica', 16, 'bold'),
            bg=self.COLORS['bg_dark'],
            fg=self.COLORS['text_primary']
        )
        title_label.pack()

        subtitle_label = tk.Label(
            header_frame,
            text="Compare Pattern Matching vs Quantitative Scoring",
            font=('Helvetica', 10),
            bg=self.COLORS['bg_dark'],
            fg=self.COLORS['text_secondary']
        )
        subtitle_label.pack()

        # Separator
        separator = tk.Frame(main_frame, height=2, bg=self.COLORS['bg_light'])
        separator.pack(fill=tk.X, pady=(10, 12))

        # ========== INPUT SECTION ==========
        input_frame = tk.Frame(main_frame, bg=self.COLORS['bg_dark'])
        input_frame.pack(fill=tk.X, pady=(0, 12))

        # Compound Name Input
        name_label = tk.Label(
            input_frame,
            text="COMPOUND NAME (Optional):",
            font=('Helvetica', 11, 'bold'),
            bg=self.COLORS['bg_dark'],
            fg=self.COLORS['text_primary'],
            anchor='w'
        )
        name_label.pack(fill=tk.X, pady=(0, 6))

        name_container = tk.Frame(
            input_frame,
            bg=self.COLORS['bg_light'],
            padx=2,
            pady=2
        )
        name_container.pack(fill=tk.X)

        self.name_entry = tk.Entry(
            name_container,
            font=('Helvetica', 11),
            bg=self.COLORS['input_bg'],
            fg=self.COLORS['input_fg'],
            insertbackground=self.COLORS['input_fg'],
            relief=tk.FLAT,
            highlightthickness=0
        )
        self.name_entry.pack(fill=tk.X, ipady=6)

        # Example name text
        name_example_label = tk.Label(
            input_frame,
            text="Example: Compound_A, Peptidomimetic_1, etc.",
            font=('Helvetica', 8),
            bg=self.COLORS['bg_dark'],
            fg=self.COLORS['text_secondary'],
            anchor='w'
        )
        name_example_label.pack(fill=tk.X, pady=(4, 0))

        # SMILES Input
        smiles_label = tk.Label(
            input_frame,
            text="PASTE SMILES DATA:",
            font=('Helvetica', 11, 'bold'),
            bg=self.COLORS['bg_dark'],
            fg=self.COLORS['text_primary'],
            anchor='w'
        )
        smiles_label.pack(fill=tk.X, pady=(12, 6))

        # Input container
        input_container = tk.Frame(
            input_frame,
            bg=self.COLORS['bg_light'],
            padx=2,
            pady=2
        )
        input_container.pack(fill=tk.X)

        self.smiles_entry = tk.Entry(
            input_container,
            font=('Courier New', 11, 'bold'),
            bg=self.COLORS['input_bg'],
            fg=self.COLORS['input_fg'],
            insertbackground=self.COLORS['input_fg'],
            relief=tk.FLAT,
            highlightthickness=0
        )
        self.smiles_entry.pack(fill=tk.X, ipady=6)

        # ========== METHOD SELECTION ==========
        method_frame = tk.Frame(
            main_frame,
            bg=self.COLORS['bg_medium'],
            padx=15,
            pady=12
        )
        method_frame.pack(fill=tk.X, pady=(0, 12))

        method_title = tk.Label(
            method_frame,
            text="SELECT PREDICTION METHOD:",
            font=('Helvetica', 10, 'bold'),
            bg=self.COLORS['bg_medium'],
            fg=self.COLORS['text_primary'],
            anchor='w'
        )
        method_title.pack(fill=tk.X, pady=(0, 8))

        # Method radio buttons
        method1 = tk.Radiobutton(
            method_frame,
            text="⚡ Simple Pattern Matching (Fast, Qualitative)",
            variable=self.method_var,
            value="pattern",
            font=('Helvetica', 9, 'bold'),
            bg=self.COLORS['bg_medium'],
            fg=self.COLORS['text_primary'],
            selectcolor=self.COLORS['bg_dark'],
            activebackground=self.COLORS['bg_medium'],
            activeforeground=self.COLORS['accent_blue']
        )
        method1.pack(anchor='w', pady=2)

        method2 = tk.Radiobutton(
            method_frame,
            text="🎯 SMARTCyp-Inspired Scoring (Advanced, Quantitative)",
            variable=self.method_var,
            value="quantitative",
            font=('Helvetica', 9, 'bold'),
            bg=self.COLORS['bg_medium'],
            fg=self.COLORS['text_primary'],
            selectcolor=self.COLORS['bg_dark'],
            activebackground=self.COLORS['bg_medium'],
            activeforeground=self.COLORS['accent_purple']
        )
        method2.pack(anchor='w', pady=2)

        method3 = tk.Radiobutton(
            method_frame,
            text="📊 Both Methods (Compare Side-by-Side)",
            variable=self.method_var,
            value="both",
            font=('Helvetica', 9, 'bold'),
            bg=self.COLORS['bg_medium'],
            fg=self.COLORS['text_primary'],
            selectcolor=self.COLORS['bg_dark'],
            activebackground=self.COLORS['bg_medium'],
            activeforeground=self.COLORS['accent_green']
        )
        method3.pack(anchor='w', pady=2)

        # ========== BUTTONS ==========
        button_frame = tk.Frame(main_frame, bg=self.COLORS['bg_dark'])
        button_frame.pack(fill=tk.X, pady=15)

        # Run button
        self.run_button = tk.Button(
            button_frame,
            text='RUN PREDICTION',
            font=('Helvetica', 12, 'bold'),
            bg=self.COLORS['btn_primary'],
            fg=self.COLORS['text_on_button'],
            activebackground=self.COLORS['btn_primary_hover'],
            activeforeground=self.COLORS['text_on_button'],
            cursor='hand2',
            relief=tk.FLAT,
            borderwidth=0,
            padx=20,
            pady=12,
            command=self.validate_and_run
        )
        self.run_button.pack(fill=tk.X, pady=(0, 8))

        # Utility buttons
        util_button_frame = tk.Frame(button_frame, bg=self.COLORS['bg_dark'])
        util_button_frame.pack(fill=tk.X)

        self.example_button = tk.Button(
            util_button_frame,
            text='LOAD EXAMPLE',
            font=('Helvetica', 9, 'bold'),
            bg=self.COLORS['btn_success'],
            fg=self.COLORS['text_on_button'],
            cursor='hand2',
            relief=tk.FLAT,
            padx=12,
            pady=8,
            command=self.load_example
        )
        self.example_button.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=(0, 4))

        self.clear_button = tk.Button(
            util_button_frame,
            text='CLEAR',
            font=('Helvetica', 9, 'bold'),
            bg=self.COLORS['btn_warning'],
            fg=self.COLORS['text_on_button'],
            cursor='hand2',
            relief=tk.FLAT,
            padx=12,
            pady=8,
            command=self.clear_input
        )
        self.clear_button.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=(4, 0))

        # ========== STATUS ==========
        status_frame = tk.Frame(
            main_frame,
            bg=self.COLORS['bg_medium'],
            padx=12,
            pady=8
        )
        status_frame.pack(fill=tk.X, pady=(0, 12))

        self.status_label = tk.Label(
            status_frame,
            text="● Ready to analyze",
            font=('Helvetica', 9, 'bold'),
            bg=self.COLORS['bg_medium'],
            fg=self.COLORS['accent_green'],
            anchor='w'
        )
        self.status_label.pack(fill=tk.X)

        # ========== INFO BOX ==========
        info_frame = tk.Frame(
            main_frame,
            bg=self.COLORS['bg_light'],
            padx=12,
            pady=10
        )
        info_frame.pack(fill=tk.X)

        info_title = tk.Label(
            info_frame,
            text="METHOD COMPARISON:",
            font=('Helvetica', 9, 'bold'),
            bg=self.COLORS['bg_light'],
            fg=self.COLORS['text_primary'],
            anchor='w'
        )
        info_title.pack(fill=tk.X, pady=(0, 5))

        info_text = [
            "• Pattern Matching: Binary detection (present/absent)",
            "• Quantitative Scoring: Ranked sites with metabolic liability scores (0-100)",
            "• Both: See differences between simple and advanced predictions"
        ]

        for line in info_text:
            lbl = tk.Label(
                info_frame,
                text=line,
                font=('Helvetica', 8),
                bg=self.COLORS['bg_light'],
                fg=self.COLORS['text_secondary'],
                anchor='w'
            )
            lbl.pack(fill=tk.X, pady=1)

    def center_window(self):
        """Center window on screen"""
        self.root.update_idletasks()
        width = self.root.winfo_width()
        height = self.root.winfo_height()
        x = (self.root.winfo_screenwidth() // 2) - (width // 2)
        y = (self.root.winfo_screenheight() // 2) - (height // 2)
        self.root.geometry(f'{width}x{height}+{x}+{y}')

    def validate_and_run(self):
        """Validate input and run analysis"""
        smiles = self.smiles_entry.get().strip()

        if not smiles:
            messagebox.showwarning(
                "Empty Input",
                "Please enter SMILES data or click 'LOAD EXAMPLE'",
                parent=self.root
            )
            return

        if len(smiles) < 2:
            messagebox.showwarning(
                "Invalid SMILES",
                "The SMILES string appears too short.",
                parent=self.root
            )
            return

        self.run_analysis_threaded(smiles)

    def run_analysis_threaded(self, smiles):
        """Run analysis in background thread"""
        thread = threading.Thread(
            target=self.run_analysis,
            args=(smiles,),
            daemon=True
        )
        thread.start()

    def run_analysis(self, smiles):
        """Execute prediction based on selected method"""
        try:
            # Update UI
            self.run_button.config(
                state='disabled',
                text='ANALYZING...',
                bg=self.COLORS['bg_medium']
            )
            self.status_label.config(
                text="● Running prediction...",
                fg=self.COLORS['accent_orange']
            )
            self.root.update()

            method = self.method_var.get()

            # Get compound name for file naming
            compound_name = self.name_entry.get().strip()
            if compound_name:
                # Sanitize filename (remove invalid characters)
                import re
                compound_name = re.sub(r'[<>:"/\\|?*]', '_', compound_name)
                pattern_filename = f"{compound_name}_pattern_result.png"
                quant_filename = f"{compound_name}_quantitative_result.png"
            else:
                pattern_filename = "pattern_matching_result.png"
                quant_filename = "quantitative_scoring_result.png"

            # Capture output
            old_stdout = sys.stdout
            sys.stdout = captured_output = io.StringIO()

            results = {}
            results['filenames'] = {'pattern': pattern_filename, 'quantitative': quant_filename}

            # Run pattern matching
            if method in ('pattern', 'both'):
                from metabolic_liability_tagger import MetabolicLiabilityTagger
                pattern_tagger = MetabolicLiabilityTagger(smiles)
                pattern_tagger.run_full_analysis(pattern_filename)
                results['pattern'] = captured_output.getvalue()

            # Run quantitative scoring
            if method in ('quantitative', 'both'):
                from smartcyp_inspired_predictor import SMARTCypInspiredPredictor
                quant_predictor = SMARTCypInspiredPredictor(smiles)
                quant_predictor.predict()
                report = quant_predictor.generate_report()
                quant_predictor.visualize(quant_filename, top_n=5)
                results['quantitative'] = report

            # Restore stdout
            sys.stdout = old_stdout

            # Update UI
            self.run_button.config(
                state='normal',
                text='RUN PREDICTION',
                bg=self.COLORS['btn_primary']
            )
            self.status_label.config(
                text="● Analysis complete!",
                fg=self.COLORS['accent_green']
            )

            # Show results
            self.show_results(results, smiles, method, compound_name)

        except Exception as e:
            self.handle_error("Analysis Error", str(e))
            import traceback
            traceback.print_exc()

    def handle_error(self, title, message):
        """Handle errors"""
        self.run_button.config(
            state='normal',
            text='RUN PREDICTION',
            bg=self.COLORS['btn_primary']
        )
        self.status_label.config(
            text=f"● Error: {title}",
            fg=self.COLORS['accent_orange']
        )
        messagebox.showerror(title, message, parent=self.root)

    def show_results(self, results, smiles, method, compound_name=""):
        """Display results in popup"""
        result_window = tk.Toplevel(self.root)

        # Set window title with compound name if provided
        if compound_name:
            result_window.title(f"Prediction Results - {compound_name}")
        else:
            result_window.title("Prediction Results")

        result_window.geometry("900x700")
        result_window.configure(bg=self.COLORS['bg_dark'])

        # Header
        header = tk.Frame(result_window, bg=self.COLORS['bg_medium'], pady=10)
        header.pack(fill=tk.X)

        # Title with compound name if provided
        title_text = f"RESULTS - {method.upper()} METHOD(S)"
        if compound_name:
            title_text += f"\nCompound: {compound_name}"

        title = tk.Label(
            header,
            text=title_text,
            font=('Helvetica', 13, 'bold'),
            bg=self.COLORS['bg_medium'],
            fg=self.COLORS['text_primary']
        )
        title.pack()

        # Results display
        text_frame = tk.Frame(result_window, bg=self.COLORS['bg_dark'], padx=10, pady=10)
        text_frame.pack(fill=tk.BOTH, expand=True)

        text_widget = scrolledtext.ScrolledText(
            text_frame,
            wrap=tk.WORD,
            font=('Courier New', 9),
            bg=self.COLORS['bg_light'],
            fg=self.COLORS['text_primary'],
            insertbackground=self.COLORS['text_primary'],
            relief=tk.FLAT
        )
        text_widget.pack(fill=tk.BOTH, expand=True)

        # Format output
        output_text = ""
        if 'pattern' in results:
            output_text += "=" * 90 + "\n"
            output_text += "PATTERN MATCHING RESULTS\n"
            output_text += "=" * 90 + "\n\n"
            output_text += results['pattern'] + "\n\n"

        if 'quantitative' in results:
            output_text += results['quantitative'] + "\n"

        text_widget.insert(tk.END, output_text)
        text_widget.config(state='disabled')

        # Buttons
        button_frame = tk.Frame(result_window, bg=self.COLORS['bg_dark'], pady=10)
        button_frame.pack(fill=tk.X)

        # View images button(s)
        filenames = results.get('filenames', {'pattern': 'pattern_matching_result.png',
                                               'quantitative': 'quantitative_scoring_result.png'})

        if method == 'both':
            view_pattern_btn = tk.Button(
                button_frame,
                text='VIEW PATTERN RESULT',
                font=('Helvetica', 10, 'bold'),
                bg=self.COLORS['btn_info'],
                fg=self.COLORS['text_on_button'],
                cursor='hand2',
                relief=tk.FLAT,
                padx=15,
                pady=8,
                command=lambda: self.view_image(filenames['pattern'])
            )
            view_pattern_btn.pack(side=tk.LEFT, padx=10, expand=True, fill=tk.X)

            view_quant_btn = tk.Button(
                button_frame,
                text='VIEW QUANTITATIVE RESULT',
                font=('Helvetica', 10, 'bold'),
                bg=self.COLORS['btn_info'],
                fg=self.COLORS['text_on_button'],
                cursor='hand2',
                relief=tk.FLAT,
                padx=15,
                pady=8,
                command=lambda: self.view_image(filenames['quantitative'])
            )
            view_quant_btn.pack(side=tk.LEFT, padx=10, expand=True, fill=tk.X)
        else:
            image_name = filenames['pattern'] if method == "pattern" else filenames['quantitative']
            view_btn = tk.Button(
                button_frame,
                text='VIEW IMAGE',
                font=('Helvetica', 10, 'bold'),
                bg=self.COLORS['btn_success'],
                fg=self.COLORS['text_on_button'],
                cursor='hand2',
                relief=tk.FLAT,
                padx=20,
                pady=10,
                command=lambda: self.view_image(image_name)
            )
            view_btn.pack(side=tk.LEFT, padx=10, expand=True, fill=tk.X)

        close_btn = tk.Button(
            button_frame,
            text='CLOSE',
            font=('Helvetica', 10, 'bold'),
            bg=self.COLORS['btn_warning'],
            fg=self.COLORS['text_on_button'],
            cursor='hand2',
            relief=tk.FLAT,
            padx=20,
            pady=10,
            command=result_window.destroy
        )
        close_btn.pack(side=tk.LEFT, padx=10, expand=True, fill=tk.X)

    def view_image(self, filename):
        """Open generated image"""
        import subprocess
        import os

        image_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename)

        if os.path.exists(image_path):
            if sys.platform == 'darwin':  # macOS
                subprocess.Popen(['open', image_path])
            elif sys.platform == 'win32':  # Windows
                os.startfile(image_path)
            else:  # Linux
                subprocess.Popen(['xdg-open', image_path])
        else:
            messagebox.showwarning(
                "Image Not Found",
                f"{filename} not found.",
                parent=self.root
            )

    def load_example(self):
        """Load example SMILES"""
        self.name_entry.delete(0, tk.END)
        self.name_entry.insert(0, "Ac-Cha-Val-Tic-Dimethyluracil")
        self.smiles_entry.delete(0, tk.END)
        self.smiles_entry.insert(0, self.default_smiles)
        self.status_label.config(
            text="● Example loaded - select method and click RUN",
            fg=self.COLORS['accent_blue']
        )

    def clear_input(self):
        """Clear input"""
        self.name_entry.delete(0, tk.END)
        self.smiles_entry.delete(0, tk.END)
        self.status_label.config(
            text="● Input cleared",
            fg=self.COLORS['text_secondary']
        )


def main():
    """Main entry point"""
    root = tk.Tk()
    app = HybridMetabolismGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
