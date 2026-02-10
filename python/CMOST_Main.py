###############################################################################
#
#     CMOST: Colon Modeling with Open Source Tool
#     created by Meher Prakash and Benjamin Misselwitz 2012 - 2016
#
#     This program is part of free software package CMOST for colo-rectal
#     cancer simulations: You can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

# POLYPCALCULATOR
#      simulation for the natural history of colorectal polyps and cancer
#      helpful for testing the effects of various screening strategies
#
#      Program created by Benjamin Misselwitz, 2010-2011

import os
import sys
import re
import copy
import pickle
import numpy as np

# Ensure the python/ directory is on the import path so sibling modules
# (calculate_sub, NumberCrunching_100000, Evaluation, etc.) can be found.
_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

import tkinter as tk
from tkinter import messagebox, filedialog, simpledialog

# ---------------------------------------------------------------------------
# INDEX CONVENTION NOTES:
#
# In the original MATLAB code, arrays are 1-based.  In this Python
# translation, all array indices are 0-based.
#
# MATLAB's  handles.Variables  is stored as  self.variables  (a dict).
# MATLAB's  handles.Variables.Screening  becomes  self.variables['Screening']
# MATLAB structs become Python dicts, consistent with NumberCrunching_100000.py,
# Evaluation.py, Colonoscopy_settings.py, Cost_Settings.py, and
# Risk_Settings.py.
#
# The MATLAB GUI was built using GUIDE (gui_mainfcn / gui_Singleton).
# This Python version uses tkinter for the GUI, following the same pattern
# as the previously converted Colonoscopy_settings.py, Cost_Settings.py,
# and Risk_Settings.py.
#
# Sub-dialogs (Colonoscopy_settings, Cost_Settings, Risk_Settings, etc.)
# are launched by importing their module-level convenience functions.
# In MATLAB, data was passed via set(0,'userdata',...)/get(0,'userdata').
# In Python, the Variables dict is passed by reference to the sub-dialog.
#
# Settings files:
#   MATLAB used .mat files loaded via load()/save().
#   Python uses .pkl (pickle) files or the auto-generated Python settings
#   modules (settings/CMOST13.py, etc.).
#
# NumberPatientsValues has 4 elements: index 0..3 (MATLAB 1..4).
# ---------------------------------------------------------------------------


# ===================================================================
#  IMPORTS OF PREVIOUSLY CONVERTED MODULES
# ===================================================================

# Import sub-dialog convenience functions from previously converted modules.
# These modules follow the pattern:  open dialog, block, return modified dict.
# If a module has not yet been converted, we provide a stub that shows a
# warning message so the rest of the GUI remains functional.

def _import_optional(module_name, func_name):
    """
    Attempt to import func_name from module_name.
    Returns the function if found, otherwise returns a stub that warns the user.
    """
    try:
        mod = __import__(module_name, fromlist=[func_name])
        return getattr(mod, func_name)
    except (ImportError, AttributeError):
        def _stub(variables, **kwargs):
            messagebox.showwarning(
                'Module Not Found',
                f'The module "{module_name}" (function "{func_name}") '
                f'has not been converted to Python yet.')
            return variables
        return _stub


# Previously converted sub-dialogs
colonoscopy_settings = _import_optional('Colonoscopy_settings', 'colonoscopy_settings')
cost_settings = _import_optional('Cost_Settings', 'cost_settings')
risk_settings = _import_optional('Risk_Settings', 'risk_settings')

# Sub-dialogs that may not yet be converted -- stubs will warn the user
screening_settings = _import_optional('Screening_Settings', 'screening_settings')
mortality_settings = _import_optional('Mortality_Settings', 'mortality_settings')
location_settings = _import_optional('Location_Settings', 'location_settings')
starter_dialog = _import_optional('Starter', 'starter')
scan_variables_dialog = _import_optional('ScanVariables', 'scan_variables')
manual_adjustments_dialog = _import_optional('ManualAdjustments', 'manual_adjustments')
sensitivity_analysis_dialog = _import_optional('Sensitivity_Analysis', 'sensitivity_analysis')

# Benchmark GUIs
step_1_benchmarks = _import_optional('Step_1_Benchmarks_EarlyAdenoma', 'step_1_benchmarks_early_adenoma')
step_2_benchmarks = _import_optional('Step_2_Benchmarks_AdvancedAdenoma', 'step_2_benchmarks_advanced_adenoma')
step_3_benchmarks = _import_optional('Step_3_Benchmarks_Carcinoma', 'step_3_benchmarks_carcinoma')
step_4_benchmarks = _import_optional('Step_4_Benchmarks_RSRCT', 'step_4_benchmarks_rsrct')

# Auto-calibration
auto_calibration_step_1 = _import_optional('Auto_Calibration_Step_1', 'auto_calibration_step_1')
auto_calibration_step_2 = _import_optional('Auto_Calibration_Step_2', 'auto_calibration_step_2')
auto_calibration_step_3 = _import_optional('Auto_Calibration_Step_3', 'auto_calibration_step_3')
auto_calibration_step_4 = _import_optional('Auto_Calibration_Step_4', 'auto_calibration_step_4')
automatic_calibration_steps123 = _import_optional('AutomaticCalibration_Steps123', 'automatic_calibration_steps123')

# Other utilities
default_benchmarks_func = _import_optional('Default_Benchmarks', 'default_benchmarks')
automatic_settings_writing_func = _import_optional('Automatic_Settings_Writing', 'automatic_settings_writing')
automatic_rs_screen_func = _import_optional('Automatic_RS_Screen', 'automatic_rs_screen')
evaluate_rs_scan_func = _import_optional('Evaluate_RS_Scan', 'evaluate_rs_scan')

# The core simulation function
try:
    from calculate_sub import calculate_sub
except ImportError:
    def calculate_sub(handles):
        messagebox.showerror('Module Not Found',
                             'calculate_sub.py not found. Cannot run simulation.')
        return handles, None


# ===================================================================
#  SAFE VARIABLE NAME CHECK  (replaces MATLAB isvarname)
# ===================================================================

def _isvarname(name):
    """
    Check whether a string is a valid variable/settings name.
    Mirrors MATLAB's isvarname: must start with a letter and contain
    only letters, digits, or underscores.
    """
    if not name:
        return False
    return bool(re.match(r'^[A-Za-z][A-Za-z0-9_]*$', name))


# ===================================================================
#  SETTINGS I/O  (replaces MATLAB load/save with .mat files)
# ===================================================================

def _load_settings(filepath):
    """
    Load a settings dictionary from a file.

    Supports:
      - .pkl (pickle) files: expected to contain a dict
      - .py (Python settings module): expected to have a 'settings' dict

    Returns
    -------
    dict or None
        The loaded settings dictionary, or None on failure.
    """
    ext = os.path.splitext(filepath)[1].lower()

    if ext == '.pkl':
        try:
            with open(filepath, 'rb') as f:
                data = pickle.load(f)
            if isinstance(data, dict):
                return data
        except Exception:
            pass
        return None

    elif ext == '.py':
        # Import the Python settings module dynamically
        try:
            import importlib.util
            spec = importlib.util.spec_from_file_location('_tmp_settings', filepath)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            if hasattr(mod, 'settings') and isinstance(mod.settings, dict):
                return copy.deepcopy(mod.settings)
        except Exception:
            pass
        return None

    elif ext == '.mat':
        # Attempt to load MATLAB .mat file using scipy if available
        try:
            import scipy.io
            mat_data = scipy.io.loadmat(filepath, squeeze_me=True)
            # The MATLAB code saves as: save(filename, 'temp')
            # so the dict key is 'temp'
            if 'temp' in mat_data:
                return _mat_struct_to_dict(mat_data['temp'])
            elif 'Variables' in mat_data:
                return _mat_struct_to_dict(mat_data['Variables'])
            else:
                # Return the first non-internal key
                for key in mat_data:
                    if not key.startswith('__'):
                        return _mat_struct_to_dict(mat_data[key])
        except ImportError:
            messagebox.showwarning('scipy not available',
                                   'Cannot load .mat files without scipy. '
                                   'Please install scipy or use .pkl/.py files.')
        except Exception:
            pass
        return None

    return None


def _mat_struct_to_dict(obj):
    """
    Recursively convert a scipy.io loaded MATLAB struct to a Python dict.
    """
    if hasattr(obj, 'dtype') and obj.dtype.names is not None:
        # structured numpy array (MATLAB struct)
        result = {}
        for name in obj.dtype.names:
            val = obj[name]
            if hasattr(val, 'item'):
                val = val.item()
            result[name] = _mat_struct_to_dict(val)
        return result
    elif isinstance(obj, np.ndarray):
        if obj.ndim == 0:
            return obj.item()
        elif obj.dtype.kind in ('U', 'S', 'O'):
            # String or object array
            if obj.size == 1:
                return str(obj.flat[0])
            return [str(x) for x in obj.flat]
        else:
            return obj.tolist()
    else:
        return obj


def _save_settings(filepath, variables):
    """
    Save a settings dictionary to a .pkl file.

    MATLAB equivalent:
        temp = handles.Variables;
        save(fullfile(pathname, filename), 'temp');

    Parameters
    ----------
    filepath : str
        Path to the output .pkl file.
    variables : dict
        The Variables dictionary to save.
    """
    with open(filepath, 'wb') as f:
        pickle.dump(variables, f, protocol=pickle.HIGHEST_PROTOCOL)


# ===================================================================
#  LOAD DEFAULT SETTINGS
# ===================================================================

def _load_default_settings(current_path):
    """
    Load the default CMOST13 settings.

    MATLAB equivalent:
        load(fullfile(CurrentPath, 'Settings', 'CMOST13.mat'))
        handles.Variables = temp;

    Tries (in order):
      1. settings/CMOST13.py (Python settings module)
      2. settings/CMOST13.pkl
      3. Settings/CMOST13.mat

    Returns
    -------
    dict
        The loaded settings dictionary.
    """
    # Try Python settings module first
    py_path = os.path.join(current_path, 'settings', 'CMOST13.py')
    if os.path.isfile(py_path):
        result = _load_settings(py_path)
        if result is not None:
            return result

    # Try pickle
    pkl_path = os.path.join(current_path, 'settings', 'CMOST13.pkl')
    if os.path.isfile(pkl_path):
        result = _load_settings(pkl_path)
        if result is not None:
            return result

    # Try .mat file in parent directory's Settings folder
    mat_path = os.path.join(os.path.dirname(current_path), 'Settings', 'CMOST13.mat')
    if os.path.isfile(mat_path):
        result = _load_settings(mat_path)
        if result is not None:
            return result

    # Fallback: try to import settings.CMOST13 as a Python module
    try:
        sys.path.insert(0, current_path)
        from settings.CMOST13 import settings
        return copy.deepcopy(settings)
    except ImportError:
        pass

    messagebox.showerror('Settings Not Found',
                         'Could not load default settings (CMOST13). '
                         'Please ensure settings/CMOST13.py exists.')
    return {}


def _load_life_table(current_path):
    """
    Load the LifeTable data.

    MATLAB equivalent:
        load(fullfile(handles.Variables.CurrentPath, 'LifeTable.mat'))
        handles.Variables.LifeTable = LifeTable;

    Tries:
      1. LifeTable.pkl in current_path
      2. LifeTable.mat in parent directory

    Returns
    -------
    numpy.ndarray or list
        The life table data.
    """
    # Try pickle first
    pkl_path = os.path.join(current_path, 'LifeTable.pkl')
    if os.path.isfile(pkl_path):
        try:
            with open(pkl_path, 'rb') as f:
                data = pickle.load(f)
            if isinstance(data, dict) and 'LifeTable' in data:
                return data['LifeTable']
            return data
        except Exception:
            pass

    # Try .mat file in parent directory
    mat_path = os.path.join(os.path.dirname(current_path), 'LifeTable.mat')
    if os.path.isfile(mat_path):
        try:
            import scipy.io
            mat_data = scipy.io.loadmat(mat_path, squeeze_me=True)
            if 'LifeTable' in mat_data:
                return mat_data['LifeTable'].tolist()
        except ImportError:
            pass
        except Exception:
            pass

    messagebox.showwarning('LifeTable Not Found',
                           'Could not load LifeTable. '
                           'Simulation may not work correctly.')
    return []


# ===================================================================
#  CMOST MAIN GUI
# ===================================================================

class CMOSTMainGUI:
    """
    Python/tkinter equivalent of the MATLAB GUIDE CMOST_Main dialog.

    This is the main entry point for the CMOST application. The GUI allows
    the user to:
      - Load/save simulation settings
      - Configure screening, surveillance, colonoscopy, cost, risk,
        location, mortality, and screening parameters via sub-dialogs
      - Set the number of patients for the simulation
      - Run single or batch simulations
      - Access benchmark calibration tools
      - Configure output options (results file, Excel, PDF)

    The Variables dictionary is the central data structure that holds all
    simulation parameters. It is passed to sub-dialogs by reference and
    modified in-place.

    MATLAB equivalent: CMOST_Main.m (GUIDE-based GUI)
    """

    def __init__(self):
        """
        Opening function.

        MATLAB equivalent: CMOST_Main_OpeningFcn
        This opening function defines variables and settings and prepares
        the graphical user interface.
        """
        # --- Build the root window -------------------------------------------
        self.root = tk.Tk()
        self.root.title('CMOST: Colon Modeling with Open Source Tool')
        self.root.configure(bg='#9999ff')  # [0.6 0.6 1] in MATLAB

        # Dictionary to hold all widget references (mirrors MATLAB handles.xxx)
        self.widgets = {}

        # %%% Name and comments
        # MATLAB:
        #   handles.Variables.Settings_Name  = 'Default';
        #   handles.Variables.Comment        = 'no comment please';
        #   handles.Variables.Identification = 'This_is_a_genuine_PolypCalculator_File';
        self.variables = {}
        self.variables['Settings_Name'] = 'Default'
        self.variables['Comment'] = 'no comment please'
        self.variables['Identification'] = 'This_is_a_genuine_PolypCalculator_File'

        # %%% the path where this program is stored
        # MATLAB:
        #   Path = mfilename('fullpath');
        #   pos = regexp(Path, [mfilename, '$']);
        #   CurrentPath = Path(1:pos-1);
        current_path = os.path.dirname(os.path.abspath(__file__))

        ########################################
        ###     Handles Variables            ###
        ########################################

        # load default settings
        # MATLAB:
        #   load(fullfile(CurrentPath, 'Settings', 'CMOST13.mat'))
        #   handles.Variables = temp;
        default_settings = _load_default_settings(current_path)
        if default_settings:
            self.variables = default_settings

        self.default_settings_filename = 'CMOST13.py'

        # MATLAB: handles.Variables.ResultsPath = fullfile(Path(1:pos-1), 'Results');
        self.variables['ResultsPath'] = os.path.join(current_path, 'Results')

        # %%% the path where this program is stored
        # MATLAB: handles.Variables.CurrentPath = CurrentPath;
        self.variables['CurrentPath'] = current_path

        # MATLAB: handles.Variables.NumberPatientsValues = [10000 25000 50000 100000];
        self.variables['NumberPatientsValues'] = [10000, 25000, 50000, 100000]

        ########################################
        ###     Life Table                  ###
        ########################################

        # we load the saved variable
        # MATLAB:
        #   load(fullfile(handles.Variables.CurrentPath, 'LifeTable.mat'))
        #   handles.Variables.LifeTable = LifeTable;
        life_table = _load_life_table(current_path)
        if life_table is not None and len(life_table) > 0:
            self.variables['LifeTable'] = life_table

        # Build the GUI
        self._build_gui()

        # Update handles structure
        # MATLAB: handles.output = hObject;
        # (not needed in Python -- the root window is self.root)

        # MATLAB: handles = MakeImagesCurrent(hObject, handles);
        self._make_images_current()

    # -----------------------------------------------------------------
    #  GUI construction
    # -----------------------------------------------------------------
    def _build_gui(self):
        """Build all widgets for the main CMOST GUI."""
        root = self.root

        # Make the window resizable and set a minimum size
        root.minsize(700, 600)

        # ===============================================================
        # Top section: Settings name, comment, default settings
        # ===============================================================
        top_frame = tk.LabelFrame(root, text='Settings', bg='#9999ff')
        top_frame.pack(padx=5, pady=5, fill='x')

        # Default settings name (display only)
        row0 = tk.Frame(top_frame, bg='#9999ff')
        row0.pack(fill='x', padx=5, pady=2)
        tk.Label(row0, text='Default Settings:', bg='#9999ff').pack(side='left')
        default_entry = tk.Entry(row0, width=30, state='disabled')
        default_entry.pack(side='left', padx=5)
        self.widgets['Default_SettingsName'] = default_entry

        # Settings name
        row1 = tk.Frame(top_frame, bg='#9999ff')
        row1.pack(fill='x', padx=5, pady=2)
        tk.Label(row1, text='Settings Name:', bg='#9999ff').pack(side='left')
        settings_name_entry = tk.Entry(row1, width=30)
        settings_name_entry.pack(side='left', padx=5)
        settings_name_entry.bind('<Return>', lambda e: self._settings_name_callback())
        self.widgets['settings_name'] = settings_name_entry

        # Comment
        row2 = tk.Frame(top_frame, bg='#9999ff')
        row2.pack(fill='x', padx=5, pady=2)
        tk.Label(row2, text='Comment:', bg='#9999ff').pack(side='left')
        comment_entry = tk.Entry(row2, width=50)
        comment_entry.pack(side='left', padx=5)
        comment_entry.bind('<Return>', lambda e: self._comment_callback())
        self.widgets['comment'] = comment_entry

        # ===============================================================
        # Save path section
        # ===============================================================
        path_frame = tk.LabelFrame(root, text='Save Path', bg='#9999ff')
        path_frame.pack(padx=5, pady=5, fill='x')

        path_row = tk.Frame(path_frame, bg='#9999ff')
        path_row.pack(fill='x', padx=5, pady=2)
        tk.Label(path_row, text='Results Path:', bg='#9999ff').pack(side='left')
        save_path_entry = tk.Entry(path_row, width=50)
        save_path_entry.pack(side='left', padx=5, fill='x', expand=True)
        save_path_entry.bind('<Return>', lambda e: self._save_data_path_callback())
        self.widgets['SaveDataPath_Edit'] = save_path_entry

        tk.Button(path_row, text='Browse...', command=self._browse_callback).pack(side='left', padx=5)

        # ===============================================================
        # Number of patients
        # ===============================================================
        patients_frame = tk.Frame(root, bg='#9999ff')
        patients_frame.pack(padx=5, pady=5, fill='x')

        tk.Label(patients_frame, text='Number of Patients:',
                 bg='#9999ff').pack(side='left')

        # MATLAB uses a popup menu (dropdown) with values [10000, 25000, 50000, 100000]
        self.num_patients_var = tk.StringVar()
        patient_options = ['10000', '25000', '50000', '100000']
        patients_menu = tk.OptionMenu(patients_frame, self.num_patients_var,
                                       *patient_options,
                                       command=lambda val: self._number_patients_callback())
        patients_menu.pack(side='left', padx=5)
        self.widgets['Number_patients'] = patients_menu

        # ===============================================================
        # Screening and Surveillance checkboxes
        # ===============================================================
        surv_frame = tk.LabelFrame(root, text='Screening and Surveillance',
                                    bg='#9999ff')
        surv_frame.pack(padx=5, pady=5, fill='x')

        # Polyp Surveillance checkbox
        self.polyp_surv_var = tk.IntVar()
        polyp_cb = tk.Checkbutton(
            surv_frame, text='Adenoma Surveillance',
            variable=self.polyp_surv_var,
            command=self._polyp_surveillance_callback,
            bg='#9999ff', activebackground='#9999ff')
        polyp_cb.pack(side='left', padx=10)
        self.widgets['Polyp_Surveillance'] = polyp_cb

        # Cancer Surveillance checkbox
        self.cancer_surv_var = tk.IntVar()
        cancer_cb = tk.Checkbutton(
            surv_frame, text='Cancer Surveillance',
            variable=self.cancer_surv_var,
            command=self._cancer_surveillance_callback,
            bg='#9999ff', activebackground='#9999ff')
        cancer_cb.pack(side='left', padx=10)
        self.widgets['Cancer_Surveillance'] = cancer_cb

        # Screening checkbox
        self.screening_var = tk.IntVar()
        screening_cb = tk.Checkbutton(
            surv_frame, text='Screening',
            variable=self.screening_var,
            command=self._screening_checkbox_callback,
            bg='#9999ff', activebackground='#9999ff')
        screening_cb.pack(side='left', padx=10)
        self.widgets['screening_checkbox'] = screening_cb

        # ===============================================================
        # Special flag
        # ===============================================================
        special_frame = tk.Frame(root, bg='#9999ff')
        special_frame.pack(padx=5, pady=2, fill='x')

        self.special_var = tk.IntVar()
        special_cb = tk.Checkbutton(
            special_frame, text='Special',
            variable=self.special_var,
            command=self._special_callback,
            bg='#9999ff', activebackground='#9999ff')
        special_cb.pack(side='left', padx=5)
        self.widgets['special'] = special_cb

        special_text_entry = tk.Entry(special_frame, width=30)
        special_text_entry.pack(side='left', padx=5)
        special_text_entry.bind('<Return>', lambda e: self._special_text_callback())
        self.widgets['special_text'] = special_text_entry

        # ===============================================================
        # Output options
        # ===============================================================
        output_frame = tk.LabelFrame(root, text='Output', bg='#9999ff')
        output_frame.pack(padx=5, pady=5, fill='x')

        # enable results
        self.results_var = tk.IntVar()
        results_cb = tk.Checkbutton(
            output_frame, text='Enable Results',
            variable=self.results_var,
            command=self._enable_results_callback,
            bg='#9999ff', activebackground='#9999ff')
        results_cb.pack(side='left', padx=10)
        self.widgets['enable_results'] = results_cb

        # enable excel
        self.excel_var = tk.IntVar()
        excel_cb = tk.Checkbutton(
            output_frame, text='Excel File',
            variable=self.excel_var,
            command=self._excel_file_callback,
            bg='#9999ff', activebackground='#9999ff')
        excel_cb.pack(side='left', padx=10)
        self.widgets['excel_file'] = excel_cb

        # enable pdf
        self.pdf_var = tk.IntVar()
        pdf_cb = tk.Checkbutton(
            output_frame, text='Enable PDF',
            variable=self.pdf_var,
            command=self._enable_pdf_callback,
            bg='#9999ff', activebackground='#9999ff')
        pdf_cb.pack(side='left', padx=10)
        self.widgets['enable_pdf'] = pdf_cb

        # ===============================================================
        # Settings buttons (sub-dialogs)
        # ===============================================================
        settings_btn_frame = tk.LabelFrame(root, text='Configure Settings',
                                            bg='#9999ff')
        settings_btn_frame.pack(padx=5, pady=5, fill='x')

        btn_row1 = tk.Frame(settings_btn_frame, bg='#9999ff')
        btn_row1.pack(fill='x', padx=5, pady=2)

        tk.Button(btn_row1, text='Colonoscopy',
                  command=self._colonoscopy_callback, width=18).pack(side='left', padx=3)
        tk.Button(btn_row1, text='Cost Settings',
                  command=self._cost_settings_callback, width=18).pack(side='left', padx=3)
        tk.Button(btn_row1, text='Risk Settings',
                  command=self._risk_settings_callback, width=18).pack(side='left', padx=3)
        tk.Button(btn_row1, text='Screening Settings',
                  command=self._screening_settings_callback, width=18).pack(side='left', padx=3)

        btn_row2 = tk.Frame(settings_btn_frame, bg='#9999ff')
        btn_row2.pack(fill='x', padx=5, pady=2)

        tk.Button(btn_row2, text='Location',
                  command=self._location_callback, width=18).pack(side='left', padx=3)
        tk.Button(btn_row2, text='Mortality Settings',
                  command=self._mortality_settings_callback, width=18).pack(side='left', padx=3)
        tk.Button(btn_row2, text='Starter',
                  command=self._starter_callback, width=18).pack(side='left', padx=3)
        tk.Button(btn_row2, text='Variable Scan',
                  command=self._variables_scan_callback, width=18).pack(side='left', padx=3)

        btn_row3 = tk.Frame(settings_btn_frame, bg='#9999ff')
        btn_row3.pack(fill='x', padx=5, pady=2)

        tk.Button(btn_row3, text='Sensitivity Analysis',
                  command=self._sensitivity_analysis_callback, width=18).pack(side='left', padx=3)
        tk.Button(btn_row3, text='Manual Adjustments',
                  command=self._manual_adjustments_callback, width=18).pack(side='left', padx=3)

        # ===============================================================
        # Benchmark buttons
        # ===============================================================
        benchmark_frame = tk.LabelFrame(root, text='Benchmark Input',
                                         bg='#9999ff')
        benchmark_frame.pack(padx=5, pady=5, fill='x')

        bm_row1 = tk.Frame(benchmark_frame, bg='#9999ff')
        bm_row1.pack(fill='x', padx=5, pady=2)

        tk.Button(bm_row1, text='Early Adenoma BM',
                  command=self._early_benchmarks_callback, width=18).pack(side='left', padx=3)
        tk.Button(bm_row1, text='Advanced Adenoma BM',
                  command=self._adv_benchmarks_callback, width=18).pack(side='left', padx=3)
        tk.Button(bm_row1, text='Carcinoma BM',
                  command=self._ca_benchmarks_callback, width=18).pack(side='left', padx=3)
        tk.Button(bm_row1, text='RSRCT BM',
                  command=self._rsrct_benchmarks_callback, width=18).pack(side='left', padx=3)

        # ===============================================================
        # Auto-calibration buttons
        # ===============================================================
        calib_frame = tk.LabelFrame(root, text='Automated Calibration',
                                     bg='#9999ff')
        calib_frame.pack(padx=5, pady=5, fill='x')

        calib_row1 = tk.Frame(calib_frame, bg='#9999ff')
        calib_row1.pack(fill='x', padx=5, pady=2)

        tk.Button(calib_row1, text='Step 1',
                  command=self._step1_callback, width=12).pack(side='left', padx=3)
        tk.Button(calib_row1, text='Step 2',
                  command=self._step2_callback, width=12).pack(side='left', padx=3)
        tk.Button(calib_row1, text='Step 3',
                  command=self._step3_callback, width=12).pack(side='left', padx=3)
        tk.Button(calib_row1, text='Step 4',
                  command=self._step4_callback, width=12).pack(side='left', padx=3)
        tk.Button(calib_row1, text='Default BM',
                  command=self._default_benchmarks_callback, width=12).pack(side='left', padx=3)

        calib_row2 = tk.Frame(calib_frame, bg='#9999ff')
        calib_row2.pack(fill='x', padx=5, pady=2)

        tk.Button(calib_row2, text='Auto Calib 1-2-3',
                  command=self._auto_calib_123_callback, width=18).pack(side='left', padx=3)
        tk.Button(calib_row2, text='Auto Calib 1-2-3 (Bootstrap)',
                  command=self._auto_calib_123_bootstrapping_callback, width=24).pack(side='left', padx=3)
        tk.Button(calib_row2, text='Auto RS',
                  command=self._automatic_rs_callback, width=12).pack(side='left', padx=3)
        tk.Button(calib_row2, text='Auto RS Reading',
                  command=self._automatic_rs_reading_callback, width=16).pack(side='left', padx=3)

        # ===============================================================
        # File I/O buttons (Load, Save)
        # ===============================================================
        io_frame = tk.Frame(root, bg='#9999ff')
        io_frame.pack(padx=5, pady=5, fill='x')

        tk.Button(io_frame, text='Load Settings',
                  command=self._load_callback, width=18).pack(side='left', padx=5)
        tk.Button(io_frame, text='Save Settings',
                  command=self._save_callback, width=18).pack(side='left', padx=5)
        tk.Button(io_frame, text='Auto Settings',
                  command=self._automatic_settings_callback, width=18).pack(side='left', padx=5)
        tk.Button(io_frame, text='Repeat Identical',
                  command=self._repeat_identical_settings_callback, width=18).pack(side='left', padx=5)

        # ===============================================================
        # Start buttons
        # ===============================================================
        start_frame = tk.Frame(root, bg='#9999ff')
        start_frame.pack(padx=5, pady=10, fill='x')

        tk.Button(start_frame, text='START', command=self._start_callback,
                  width=20, height=2, bg='green', fg='white',
                  font=('TkDefaultFont', 12, 'bold')).pack(side='left', padx=10)

        start_batch_btn = tk.Button(start_frame, text='START BATCH',
                                     command=self._start_batch_callback,
                                     width=20, height=2)
        start_batch_btn.pack(side='left', padx=10)
        self.widgets['start_batch'] = start_batch_btn

    ########################################
    ###     Make Images Current          ###
    ########################################

    def _make_images_current(self):
        """
        This function is called whenever a change is made within the GUI.
        This function makes these changes visible.

        MATLAB equivalent: MakeImagesCurrent(hObject, handles)
        """
        variables = self.variables

        # MATLAB: tmp = find(handles.Variables.NumberPatientsValues == handles.Variables.Number_patients);
        #         set(handles.Number_patients, 'value', tmp)
        num_patients = variables.get('Number_patients', 100000)
        num_patients_values = variables.get('NumberPatientsValues', [10000, 25000, 50000, 100000])
        if num_patients in num_patients_values:
            self.num_patients_var.set(str(num_patients))
        else:
            self.num_patients_var.set(str(num_patients_values[-1]))

        # MATLAB: set(handles.SaveDataPath_Edit, 'string', handles.Variables.ResultsPath)
        self._set_entry('SaveDataPath_Edit', variables.get('ResultsPath', ''))

        # MATLAB: set(handles.settings_name, 'string', handles.Variables.Settings_Name)
        self._set_entry('settings_name', variables.get('Settings_Name', ''))

        # MATLAB: set(handles.comment, 'string', handles.Variables.Comment)
        self._set_entry('comment', variables.get('Comment', ''))

        # MATLAB: if isequal(handles.Variables.Polyp_Surveillance, 'off') ...
        if variables.get('Polyp_Surveillance', 'off') == 'off':
            self.polyp_surv_var.set(0)
        else:
            self.polyp_surv_var.set(1)

        # MATLAB: if isequal(handles.Variables.Cancer_Surveillance, 'off') ...
        if variables.get('Cancer_Surveillance', 'off') == 'off':
            self.cancer_surv_var.set(0)
        else:
            self.cancer_surv_var.set(1)

        # MATLAB: if isequal(handles.Variables.Screening.Mode, 'off') ...
        screening = variables.get('Screening', {})
        if screening.get('Mode', 'off') == 'off':
            self.screening_var.set(0)
        else:
            self.screening_var.set(1)

        # MATLAB: set(handles.SaveDataPath_Edit, 'string', ...)  (duplicate in original)
        # MATLAB: set(handles.settings_name, 'string', ...)  (duplicate in original)
        # MATLAB: set(handles.comment, 'string', ...)  (duplicate in original)
        # (already set above)

        # MATLAB: if isequal(handles.Variables.Starter.CurrentSummary, 'none')
        #             set(handles.start_batch, 'enable', 'off')
        #         else
        #             set(handles.start_batch, 'enable', 'on')
        #         end
        starter = variables.get('Starter', {})
        current_summary = starter.get('CurrentSummary', 'none')
        if current_summary == 'none' or current_summary is None:
            self.widgets['start_batch'].config(state='disabled')
        else:
            self.widgets['start_batch'].config(state='normal')

        # MATLAB: if isequal(handles.Variables.SpecialFlag, 'off') ...
        if variables.get('SpecialFlag', 'off') == 'off':
            self.special_var.set(0)
            self.widgets['special_text'].config(state='disabled')
        else:
            self.special_var.set(1)
            self.widgets['special_text'].config(state='normal')

        # MATLAB: if handles.Variables.ResultsFlag ...
        if variables.get('ResultsFlag', False):
            self.results_var.set(1)
        else:
            self.results_var.set(0)

        # MATLAB: if handles.Variables.ExcelFlag ...
        if variables.get('ExcelFlag', False):
            self.excel_var.set(1)
        else:
            self.excel_var.set(0)

        # MATLAB: if handles.Variables.DispFlag ...
        if variables.get('DispFlag', False):
            self.pdf_var.set(1)
        else:
            self.pdf_var.set(0)

        # MATLAB: set(handles.special_text, 'string', handles.Variables.SpecialText)
        self._set_entry('special_text', variables.get('SpecialText', ''))

        # MATLAB: set(handles.Default_SettingsName, 'string', filename, 'enable', 'off')
        default_entry = self.widgets['Default_SettingsName']
        default_entry.config(state='normal')
        default_entry.delete(0, tk.END)
        default_entry.insert(0, self.default_settings_filename)
        default_entry.config(state='disabled')

    # -----------------------------------------------------------------
    #  Helper methods for widget manipulation
    # -----------------------------------------------------------------

    def _set_entry(self, widget_name, text):
        """Helper: set the text of an Entry widget."""
        entry = self.widgets[widget_name]
        state = entry.cget('state')
        if state == 'disabled':
            entry.config(state='normal')
            entry.delete(0, tk.END)
            entry.insert(0, str(text))
            entry.config(state='disabled')
        else:
            entry.delete(0, tk.END)
            entry.insert(0, str(text))

    def _get_entry(self, widget_name):
        """Helper: get the text content of an Entry widget."""
        return self.widgets[widget_name].get()

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%%       Screening and Surveillance                  %%%
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # enable adenoma surveillance
    def _polyp_surveillance_callback(self):
        """
        MATLAB equivalent: Polyp_Surveillance_Callback

        tmp = get(handles.Polyp_Surveillance, 'value');
        if isequal(tmp, 1)
            handles.Variables.Polyp_Surveillance = 'on';
        else
            handles.Variables.Polyp_Surveillance = 'off';
        end
        handles = MakeImagesCurrent(hObject, handles);
        """
        tmp = self.polyp_surv_var.get()
        if tmp == 1:
            self.variables['Polyp_Surveillance'] = 'on'
        else:
            self.variables['Polyp_Surveillance'] = 'off'
        self._make_images_current()

    # enable cancer surveillance
    def _cancer_surveillance_callback(self):
        """
        MATLAB equivalent: Cancer_Surveillance_Callback

        tmp = get(handles.Cancer_Surveillance, 'value');
        if isequal(tmp, 1)
            handles.Variables.Cancer_Surveillance = 'on';
        else
            handles.Variables.Cancer_Surveillance = 'off';
        end
        handles = MakeImagesCurrent(hObject, handles);
        """
        tmp = self.cancer_surv_var.get()
        if tmp == 1:
            self.variables['Cancer_Surveillance'] = 'on'
        else:
            self.variables['Cancer_Surveillance'] = 'off'
        self._make_images_current()

    # enable screening
    def _screening_checkbox_callback(self):
        """
        MATLAB equivalent: screening_checkbox_Callback

        tmp = get(handles.screening_checkbox, 'value');
        if isequal(tmp, 1)
            handles.Variables.Screening.Mode = 'on';
        else
            handles.Variables.Screening.Mode = 'off';
        end
        handles = MakeImagesCurrent(hObject, handles);
        """
        tmp = self.screening_var.get()
        if 'Screening' not in self.variables:
            self.variables['Screening'] = {}
        if tmp == 1:
            self.variables['Screening']['Mode'] = 'on'
        else:
            self.variables['Screening']['Mode'] = 'off'
        self._make_images_current()

    # Special
    def _special_callback(self):
        """
        MATLAB equivalent: special_Callback

        tmp = get(handles.special, 'value');
        if isequal(tmp, 1)
            handles.Variables.SpecialFlag = 'on';
        else
            handles.Variables.SpecialFlag = 'off';
        end
        handles = MakeImagesCurrent(hObject, handles);
        """
        tmp = self.special_var.get()
        if tmp == 1:
            self.variables['SpecialFlag'] = 'on'
        else:
            self.variables['SpecialFlag'] = 'off'
        self._make_images_current()

    # special text
    def _special_text_callback(self):
        """
        MATLAB equivalent: special_text_Callback

        handles.Variables.SpecialText = get(handles.special_text, 'string');
        handles = MakeImagesCurrent(hObject, handles);
        """
        self.variables['SpecialText'] = self._get_entry('special_text')
        self._make_images_current()

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%%         OUTPUT                                    %%%
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # enable results
    def _enable_results_callback(self):
        """
        MATLAB equivalent: enable_results_Callback

        tmp = get(handles.enable_results, 'value');
        if isequal(tmp, 1)
            handles.Variables.ResultsFlag = true;
        else
            handles.Variables.ResultsFlag = false;
        end
        handles = MakeImagesCurrent(hObject, handles);
        """
        tmp = self.results_var.get()
        if tmp == 1:
            self.variables['ResultsFlag'] = True
        else:
            self.variables['ResultsFlag'] = False
        self._make_images_current()

    # enable excel
    def _excel_file_callback(self):
        """
        MATLAB equivalent: excel_file_Callback

        tmp = get(handles.excel_file, 'value');
        if isequal(tmp, 1)
            handles.Variables.ExcelFlag = true;
        else
            handles.Variables.ExcelFlag = false;
        end
        handles = MakeImagesCurrent(hObject, handles);
        """
        tmp = self.excel_var.get()
        if tmp == 1:
            self.variables['ExcelFlag'] = True
        else:
            self.variables['ExcelFlag'] = False
        self._make_images_current()

    # enable pdf
    def _enable_pdf_callback(self):
        """
        MATLAB equivalent: enable_pdf_Callback

        tmp = get(handles.enable_pdf, 'value');
        if isequal(tmp, 1)
            handles.Variables.DispFlag = true;
        else
            handles.Variables.DispFlag = false;
        end
        handles = MakeImagesCurrent(hObject, handles);
        """
        tmp = self.pdf_var.get()
        if tmp == 1:
            self.variables['DispFlag'] = True
        else:
            self.variables['DispFlag'] = False
        self._make_images_current()

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%%         Paths and settings                        %%%
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # settings
    def _settings_name_callback(self):
        """
        MATLAB equivalent: settings_name_Callback

        tmp = get(handles.settings_name, 'string');
        if isvarname(tmp)
            handles.Variables.Settings_Name = tmp;
            handles = MakeImagesCurrent(hObject, handles);
        else
            errordlg('Invalid name')
        end
        """
        tmp = self._get_entry('settings_name')
        if _isvarname(tmp):
            self.variables['Settings_Name'] = tmp
            self._make_images_current()
        else:
            messagebox.showerror('Error', 'Invalid name')

    # save data path
    def _save_data_path_callback(self):
        """
        MATLAB equivalent: SaveDataPath_Edit_Callback

        tmp = get(handles.SaveDataPath_Edit, 'string');
        if isdir(tmp)
            handles.Variables.ResultsPath = tmp;
            handles = MakeImagesCurrent(hObject, handles);
        else
            errordlg('Invalid name')
        end
        """
        tmp = self._get_entry('SaveDataPath_Edit')
        if os.path.isdir(tmp):
            self.variables['ResultsPath'] = tmp
            self._make_images_current()
        else:
            messagebox.showerror('Error', 'Invalid name')

    # browse
    def _browse_callback(self):
        """
        MATLAB equivalent: Browse_Callback

        folder_name = uigetdir(handles.Variables.ResultsPath, 'Select folder to save files');
        if ~isequal(folder_name, 0)
            handles.Variables.ResultsPath = folder_name;
            guidata(hObject, handles)
            handles = MakeImagesCurrent(hObject, handles);
        end
        """
        folder_name = filedialog.askdirectory(
            initialdir=self.variables.get('ResultsPath', ''),
            title='Select folder to save files')
        if folder_name:
            self.variables['ResultsPath'] = folder_name
            self._make_images_current()

    # comment
    def _comment_callback(self):
        """
        MATLAB equivalent: comment_Callback

        tmp = get(handles.comment, 'string');
        handles.Variables.Comment = tmp;
        handles = MakeImagesCurrent(hObject, handles);
        """
        tmp = self._get_entry('comment')
        self.variables['Comment'] = tmp
        self._make_images_current()

    # --- Executes on button press in Load.
    def _load_callback(self):
        """
        MATLAB equivalent: Load_Callback

        [filename, pathname] = uigetfile('.mat', 'Loading cell data', handles.Variables.CurrentPath);
        if isequal(filename, 0), return, end
        temp = importdata(fullfile(pathname, filename));
        ...  (validation and loading logic)
        """
        filepath = filedialog.askopenfilename(
            initialdir=self.variables.get('CurrentPath', ''),
            title='Loading cell data',
            filetypes=[('Settings files', '*.pkl *.py *.mat'),
                       ('Pickle files', '*.pkl'),
                       ('Python files', '*.py'),
                       ('MATLAB files', '*.mat'),
                       ('All files', '*.*')])
        if not filepath:
            return

        filename = os.path.basename(filepath)
        pathname = os.path.dirname(filepath)

        temp = _load_settings(filepath)

        if isinstance(temp, dict):
            if temp.get('Identification') == 'This_is_a_genuine_PolypCalculator_File':
                # Valid file
                self.variables = temp
                self.variables['NumberPatientsValues'] = [10000, 25000, 50000, 100000]

                # MATLAB: if ~isfield(handles.Variables, 'MortalityCorrection')
                #             handles.Variables.MortalityCorrectionGraph(1:150) = 1;
                #         end
                if 'MortalityCorrection' not in self.variables:
                    self.variables['MortalityCorrectionGraph'] = [1.0] * 150

                # MATLAB: if ~isfield(handles.Variables,'DwellSpeed')
                #             handles.Variables.DwellSpeed = 'Fast';
                #         end
                if 'DwellSpeed' not in self.variables:
                    self.variables['DwellSpeed'] = 'Fast'

                flag = 0
                try:
                    self._make_images_current()
                    # MATLAB: if ~isfield(handles.Variables, 'RiskCorrelation')
                    #             handles.Variables.RiskCorrelation = 'on';
                    #         end
                    if 'RiskCorrelation' not in self.variables:
                        self.variables['RiskCorrelation'] = 'on'
                except Exception:
                    flag = 1
            else:
                flag = 1
        else:
            flag = 1

        if flag == 1:
            messagebox.showerror('Error', 'Loaded settings not valid')
        else:
            messagebox.showinfo('Success', 'Settings successfully loaded')
            self.variables['CurrentPath'] = pathname
            self.default_settings_filename = filename
            self._make_images_current()

    # save
    def _save_callback(self):
        """
        MATLAB equivalent: Save_Callback

        [filename, pathname] = uiputfile('.mat', 'Saving settings', ...
            fullfile(handles.Variables.CurrentPath, [handles.Variables.Settings_Name, '.mat']));
        ... (saving logic)
        """
        initial_filename = self.variables.get('Settings_Name', 'Default') + '.pkl'
        filepath = filedialog.asksaveasfilename(
            initialdir=self.variables.get('CurrentPath', ''),
            initialfile=initial_filename,
            title='Saving settings',
            filetypes=[('Pickle files', '*.pkl'),
                       ('All files', '*.*')],
            defaultextension='.pkl')
        if not filepath:
            return

        filename = os.path.basename(filepath)
        pathname = os.path.dirname(filepath)

        # MATLAB: if ischar(filename)
        #             if isempty(regexp(filename, '.mat$', 'once'))
        #                 filename = strcat(filename, '.mat');
        #             end
        if not filename.endswith('.pkl'):
            filename = filename + '.pkl'
            filepath = os.path.join(pathname, filename)

        try:
            # MATLAB: handles.Variables.Settings_Name = strrep(filename, '.mat', '');
            self.variables['Settings_Name'] = filename.replace('.pkl', '')

            _save_settings(filepath, self.variables)

            if not os.path.isfile(filepath):
                messagebox.showerror('Error', 'Saving settings failed')
            else:
                messagebox.showinfo('Success', 'Settings successfully saved')
                self.variables['CurrentPath'] = pathname
                self._make_images_current()
        except Exception:
            messagebox.showerror(
                'Error',
                'Problems saving settings. You might not have permission to '
                'write in this place or the contact to the server might be '
                'lost. Please check connections and try again.')

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%%  Sub-dialog callbacks                             %%%
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def _starter_callback(self):
        """
        MATLAB equivalent: Starter_Callback

        set(0, 'userdata', handles.Variables);
        uiwait(Starter)
        handles.Variables = get(0, 'userdata');
        handles = MakeImagesCurrent(hObject, handles);
        """
        self.variables = starter_dialog(self.variables)
        self._make_images_current()

    def _variables_scan_callback(self):
        """
        MATLAB equivalent: VariablesScan_Callback

        set(0, 'userdata', handles.Variables);
        uiwait(ScanVariables)
        handles.Variables = get(0, 'userdata');
        handles = MakeImagesCurrent(hObject, handles);
        """
        self.variables = scan_variables_dialog(self.variables)
        self._make_images_current()

    def _colonoscopy_callback(self):
        """
        MATLAB equivalent: Colonoscopy_Callback

        set(0, 'userdata', handles.Variables);
        uiwait(Colonoscopy_settings)
        handles.Variables = get(0, 'userdata');
        handles = MakeImagesCurrent(hObject, handles);
        """
        colonoscopy_settings(self.variables)
        self._make_images_current()

    def _location_callback(self):
        """
        MATLAB equivalent: Location_Callback

        set(0, 'userdata', handles.Variables);
        uiwait(Location_Settings)
        handles.Variables = get(0, 'userdata');
        handles = MakeImagesCurrent(hObject, handles);
        """
        self.variables = location_settings(self.variables)
        self._make_images_current()

    def _cost_settings_callback(self):
        """
        MATLAB equivalent: Cost_settings_Callback

        set(0, 'userdata', handles.Variables);
        uiwait(Cost_Settings)
        handles.Variables = get(0, 'userdata');
        handles = MakeImagesCurrent(hObject, handles);
        """
        cost_settings(self.variables)
        self._make_images_current()

    def _risk_settings_callback(self):
        """
        MATLAB equivalent: Risk_Settings_Callback

        set(0, 'userdata', handles.Variables);
        uiwait(Risk_Settings)
        handles.Variables = get(0, 'userdata');
        handles = MakeImagesCurrent(hObject, handles);
        """
        risk_settings(self.variables)
        self._make_images_current()

    def _screening_settings_callback(self):
        """
        MATLAB equivalent: Screening_Settings_Callback

        set(0, 'userdata', handles.Variables);
        uiwait(Screening_Settings)
        handles.Variables = get(0, 'userdata');
        handles = MakeImagesCurrent(hObject, handles);
        """
        self.variables = screening_settings(self.variables)
        self._make_images_current()

    def _mortality_settings_callback(self):
        """
        MATLAB equivalent: Mortality_Settings_Callback

        set(0, 'userdata', handles.Variables);
        uiwait(Mortality_Settings)
        handles.Variables = get(0, 'userdata');
        handles = MakeImagesCurrent(hObject, handles);
        """
        self.variables = mortality_settings(self.variables)
        self._make_images_current()

    def _sensitivity_analysis_callback(self):
        """
        MATLAB equivalent: Sensitivity_Analysis_Callback

        set(0, 'userdata', handles.Variables);
        Sensitivity_Analysis
        handles.Variables = get(0, 'userdata');
        handles = MakeImagesCurrent(hObject, handles);
        """
        self.variables = sensitivity_analysis_dialog(self.variables)
        self._make_images_current()

    def _manual_adjustments_callback(self):
        """
        MATLAB equivalent: ManualAdjustments_Callback

        set(0, 'userdata', handles.Variables);
        uiwait(ManualAdjustments)
        handles.Variables = get(0, 'userdata');
        handles = MakeImagesCurrent(hObject, handles);
        """
        self.variables = manual_adjustments_dialog(self.variables)
        self._make_images_current()

    def _number_patients_callback(self):
        """
        MATLAB equivalent: Number_patients_Callback

        tmp = get(handles.Number_patients, 'value');
        value = handles.Variables.NumberPatientsValues(tmp);
        handles.Variables.Number_patients = value;
        handles = MakeImagesCurrent(hObject, handles);
        """
        try:
            value = int(self.num_patients_var.get())
            self.variables['Number_patients'] = value
        except ValueError:
            pass
        self._make_images_current()

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #   functions for benchmark input                          %
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # opens GUI to adjust early adenoma benchmarks
    def _early_benchmarks_callback(self):
        """
        MATLAB equivalent: EarlyBenchmarks_Callback

        set(0, 'userdata', handles.Variables);
        Step_1_Benchmarks_EarlyAdenoma
        handles.Variables = get(0, 'userdata');
        handles = MakeImagesCurrent(hObject, handles);
        """
        self.variables = step_1_benchmarks(self.variables)
        self._make_images_current()

    # opens GUI to adjust advanced adenoma benchmarks including dwell time
    def _adv_benchmarks_callback(self):
        """
        MATLAB equivalent: AdvBenchmarks_Callback

        set(0, 'userdata', handles.Variables);
        Step_2_Benchmarks_AdvancedAdenoma
        handles.Variables = get(0, 'userdata');
        handles = MakeImagesCurrent(hObject, handles);
        """
        self.variables = step_2_benchmarks(self.variables)
        self._make_images_current()

    # opens GUI to adjust carcinoma benchmarks including polyp danger
    def _ca_benchmarks_callback(self):
        """
        MATLAB equivalent: CaBenchmarks_Callback

        set(0, 'userdata', handles.Variables);
        Step_3_Benchmarks_Carcinoma
        handles.Variables = get(0, 'userdata');
        handles = MakeImagesCurrent(hObject, handles);
        """
        self.variables = step_3_benchmarks(self.variables)
        self._make_images_current()

    # opens GUI to adjust benchmarks for direct cancer, using the
    # rectosigmoidoscopy study for benchmarking
    def _rsrct_benchmarks_callback(self):
        """
        MATLAB equivalent: RSRCT_benchmarks_Callback

        set(0, 'userdata', handles.Variables);
        Step_4_Benchmarks_RSRCT
        handles.Variables = get(0, 'userdata');
        handles = MakeImagesCurrent(hObject, handles);
        """
        self.variables = step_4_benchmarks(self.variables)
        self._make_images_current()

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #   functions for automated adjustment of program settings %
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # --- Executes on button press in Step1.
    def _step1_callback(self):
        """
        MATLAB equivalent: Step1_Callback

        set(0, 'userdata', handles.Variables);
        Auto_Calibration_Step_1
        handles.Variables = get(0, 'userdata');
        handles = MakeImagesCurrent(hObject, handles);
        """
        self.variables = auto_calibration_step_1(self.variables)
        self._make_images_current()

    def _step2_callback(self):
        """
        MATLAB equivalent: Step2_Callback

        set(0, 'userdata', handles.Variables);
        Auto_Calibration_Step_2
        handles.Variables = get(0, 'userdata');
        handles = MakeImagesCurrent(hObject, handles);
        """
        self.variables = auto_calibration_step_2(self.variables)
        self._make_images_current()

    def _step3_callback(self):
        """
        MATLAB equivalent: Step3_Callback

        set(0, 'userdata', handles.Variables);
        Auto_Calibration_Step_3
        handles.Variables = get(0, 'userdata');
        handles = MakeImagesCurrent(hObject, handles);
        """
        self.variables = auto_calibration_step_3(self.variables)
        self._make_images_current()

    def _step4_callback(self):
        """
        MATLAB equivalent: Step4_Callback

        set(0, 'userdata', handles.Variables);
        Auto_Calibration_Step_4
        handles.Variables = get(0, 'userdata');
        handles = MakeImagesCurrent(hObject, handles);
        """
        self.variables = auto_calibration_step_4(self.variables)
        self._make_images_current()

    # Default Benchmarks.
    def _default_benchmarks_callback(self):
        """
        MATLAB equivalent: DefaultBenchmarks_Callback

        handles = Default_Benchmarks(handles);
        msgbox('Default settings restored')
        """
        # MATLAB: handles = Default_Benchmarks(handles);
        # The function takes the full handles struct; we pass variables
        self.variables = default_benchmarks_func(self.variables)
        messagebox.showinfo('Info', 'Default settings restored')

    # Automatic settings; routine writes a series of files which test basic
    # characteristics of a given set of settings
    def _automatic_settings_callback(self):
        """
        MATLAB equivalent: Automatic_settings_Callback

        Automatic_Settings_Writing(handles)
        """
        # MATLAB passes the full handles struct; we wrap in a dict
        handles_wrapper = {'Variables': self.variables}
        automatic_settings_writing_func(handles_wrapper)

    def _repeat_identical_settings_callback(self):
        """
        MATLAB equivalent: Repeate_identical_settings_Callback

        PathName = uigetdir(handles.Variables.ResultsPath, 'Select folder to save files');
        if isequal(PathName, 0), return, end

        BaseName = handles.Variables.Settings_Name;
        answer = inputdlg({'Please select base name for files!', 'Number replicas'},
                          'Write settings repetitively', 1, {BaseName, num2str(10)});
        if isempty(answer), return, end
        BaseName = answer{1};
        [NumberReplica, status] = str2num(answer{2});
        if ~status, return, end

        Variables = handles.Variables;
        for f = 1:NumberReplica
            ... (naming and saving logic)
        end
        msgbox('files have been written')
        """
        path_name = filedialog.askdirectory(
            initialdir=self.variables.get('ResultsPath', ''),
            title='Select folder to save files')
        if not path_name:
            return

        base_name = self.variables.get('Settings_Name', 'Default')

        # Python equivalent of inputdlg: use simpledialog
        base_name_input = simpledialog.askstring(
            'Write settings repetitively',
            'Please select base name for files:',
            initialvalue=base_name)
        if not base_name_input:
            return

        number_replica_str = simpledialog.askstring(
            'Write settings repetitively',
            'Number of replicas:',
            initialvalue='10')
        if not number_replica_str:
            return

        try:
            number_replica = int(number_replica_str)
        except ValueError:
            return

        # MATLAB naming logic:
        # if and(NumberReplica>100, f<10)
        #     tmp = [BaseName '_00' num2str(f)];
        # elseif and(NumberReplica>100, f<100)
        #     tmp = [BaseName '_0' num2str(f)];
        # elseif f<10
        #     tmp = [BaseName '_0' num2str(f)];
        # else
        #     tmp = [BaseName '_' num2str(f)];
        # end
        for f in range(1, number_replica + 1):  # MATLAB f=1:NumberReplica
            if number_replica > 100 and f < 10:
                tmp = f'{base_name_input}_00{f}'
            elif number_replica > 100 and f < 100:
                tmp = f'{base_name_input}_0{f}'
            elif f < 10:
                tmp = f'{base_name_input}_0{f}'
            else:
                tmp = f'{base_name_input}_{f}'

            _save_settings(os.path.join(path_name, tmp + '.pkl'), self.variables)

        messagebox.showinfo('Info', 'Files have been written')

    def _automatic_rs_callback(self):
        """
        MATLAB equivalent: Automatic_RS_Callback

        Automatic_RS_Screen(handles)
        """
        handles_wrapper = {'Variables': self.variables}
        automatic_rs_screen_func(handles_wrapper)

    def _automatic_rs_reading_callback(self):
        """
        MATLAB equivalent: Automatic_RS_reading_Callback

        [AnalysisPipeline, DirectCancerSpeed, mod] = Evaluate_RS_Scan(handles.Variables.ResultsPath);
        if isequal(mod, 'OK')
            button = questdlg(sprintf('Direct cancer speed %.2f. Do you want to keep?', ...
                DirectCancerSpeed), 'title');
            if isequal(button, 'Yes')
                handles.Variables.DirectCancerSpeed = DirectCancerSpeed/10000000;
                handles.Variables.ResultsPath = AnalysisPipeline;
                handles = MakeImagesCurrent(hObject, handles);
            end
        end
        """
        try:
            analysis_pipeline, direct_cancer_speed, mod_result = evaluate_rs_scan_func(
                self.variables.get('ResultsPath', ''))
        except Exception:
            messagebox.showerror('Error', 'Failed to evaluate RS scan.')
            return

        if mod_result == 'OK':
            answer = messagebox.askyesno(
                'title',
                f'Direct cancer speed {direct_cancer_speed:.2f}. '
                f'Do you want to keep the calculated value for direct cancer?')
            if answer:
                self.variables['DirectCancerSpeed'] = direct_cancer_speed / 10000000
                self.variables['ResultsPath'] = analysis_pipeline
                self._make_images_current()

    def _auto_calib_123_callback(self):
        """
        MATLAB equivalent: Auto_calib_123_Callback

        handles = AutomaticCalibration_Steps123(handles, 'normal');
        handles = MakeImagesCurrent(hObject, handles);
        """
        handles_wrapper = {'Variables': self.variables}
        handles_wrapper = automatic_calibration_steps123(handles_wrapper, 'normal')
        if isinstance(handles_wrapper, dict) and 'Variables' in handles_wrapper:
            self.variables = handles_wrapper['Variables']
        self._make_images_current()

    def _auto_calib_123_bootstrapping_callback(self):
        """
        MATLAB equivalent: Auto_calib_123_bootstrapping_Callback

        handles = AutomaticCalibration_Steps123(handles, 'bootstrapping');
        handles = MakeImagesCurrent(hObject, handles);
        """
        handles_wrapper = {'Variables': self.variables}
        handles_wrapper = automatic_calibration_steps123(handles_wrapper, 'bootstrapping')
        if isinstance(handles_wrapper, dict) and 'Variables' in handles_wrapper:
            self.variables = handles_wrapper['Variables']
        self._make_images_current()

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #   START                                               %
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def _start_callback(self):
        """
        MATLAB equivalent: Start_Callback

        [handles, BM] = CalculateSub(handles);
        """
        handles_wrapper = {'Variables': self.variables}
        handles_wrapper, bm = calculate_sub(handles_wrapper)
        if isinstance(handles_wrapper, dict) and 'Variables' in handles_wrapper:
            self.variables = handles_wrapper['Variables']
        self._make_images_current()

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%%         START BATCH                               %%%
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def _start_batch_callback(self):
        """
        MATLAB equivalent: start_batch_Callback

        OldVariables = handles.Variables;
        if isdir(handles.Variables.ResultsPath)
            ResultsPath = handles.Variables.ResultsPath;
        elseif isdir(handles.Variables.CurrentPath)
            ResultsPath = handles.Variables.CurrentPath;
        else
            ResultsPath = uigetdir(handles.Variables.CurrentPath, 'Please select path to save results');
        end
        handles.Variables.StarterFlag = 'on';

        try
            for f = 1 : length(handles.Variables.Starter.CurrentSummary)
                handles.Variables = importdata(fullfile(...));
                handles.Variables.StarterFlag = 'on';
                handles.Variables.Starter = OldVariables.Starter;
                handles.Variables.Starter.Counter = f;
                handles.Variables.ResultsPath = ResultsPath;
                [handles, BM] = CalculateSub(handles);
            end
            handles.Variables = OldVariables;
        catch
            handles.Variables = OldVariables;
            handles = MakeImagesCurrent(hObject, handles);
            rethrow(lasterror)
        end
        """
        old_variables = copy.deepcopy(self.variables)

        if os.path.isdir(self.variables.get('ResultsPath', '')):
            results_path = self.variables['ResultsPath']
        elif os.path.isdir(self.variables.get('CurrentPath', '')):
            results_path = self.variables['CurrentPath']
        else:
            results_path = filedialog.askdirectory(
                initialdir=self.variables.get('CurrentPath', ''),
                title='Please select path to save results')
            if not results_path:
                return

        self.variables['StarterFlag'] = 'on'

        try:
            starter = self.variables.get('Starter', {})
            current_summary = starter.get('CurrentSummary', [])
            current_path_list = starter.get('CurrentPath', [])

            # MATLAB: for f = 1 : length(handles.Variables.Starter.CurrentSummary)
            if isinstance(current_summary, list):
                num_files = len(current_summary)
            else:
                num_files = 0

            for f in range(num_files):  # 0-based (MATLAB f=1:length)
                # MATLAB: handles.Variables = importdata(fullfile(...CurrentPath{f}, ...CurrentSummary{f}));
                file_path = os.path.join(current_path_list[f], current_summary[f])
                loaded = _load_settings(file_path)
                if loaded is not None:
                    self.variables = loaded
                else:
                    raise RuntimeError(f'Could not load settings from: {file_path}')

                self.variables['StarterFlag'] = 'on'
                self.variables['Starter'] = copy.deepcopy(old_variables['Starter'])
                # MATLAB: handles.Variables.Starter.Counter = f;
                # (1-based in MATLAB)
                self.variables['Starter']['Counter'] = f + 1
                self.variables['ResultsPath'] = results_path

                # MATLAB: [handles, BM] = CalculateSub(handles);
                handles_wrapper = {'Variables': self.variables}
                handles_wrapper, bm = calculate_sub(handles_wrapper)
                if isinstance(handles_wrapper, dict) and 'Variables' in handles_wrapper:
                    self.variables = handles_wrapper['Variables']

            self.variables = old_variables

        except Exception as e:
            self.variables = old_variables
            self._make_images_current()
            messagebox.showerror('Error', str(e))
            raise

    # -----------------------------------------------------------------
    #  Run  (blocking main loop)
    # -----------------------------------------------------------------
    def run(self):
        """Start the tkinter main loop (blocking)."""
        self.root.mainloop()


# ===================================================================
#  MODULE-LEVEL CONVENIENCE FUNCTION
# ===================================================================

def cmost_main():
    """
    Open the CMOST Main GUI.

    This is the main entry point for the application, equivalent to
    calling CMOST_Main in MATLAB which opens the GUIDE-based GUI.
    """
    gui = CMOSTMainGUI()
    gui.run()


# ===================================================================
#  STANDALONE EXECUTION
# ===================================================================

if __name__ == '__main__':
    cmost_main()
