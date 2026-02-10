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

import numpy as np
import copy
import tkinter as tk
from tkinter import messagebox

# ---------------------------------------------------------------------------
# INDEX CONVENTION NOTES:
#
# In the original MATLAB code, arrays are 1-based.  In this Python
# translation, all array indices are 0-based.
#
# ColonoscopyRate has 20 elements: index 0..19 (MATLAB 1..20).
# ColonoscopyLikelyhood has 150 elements: index 0..149 (MATLAB 1..150).
# Colo_Detection has 10 elements: index 0..9 (MATLAB 1..10).
# RectoSigmo_Detection has 10 elements: index 0..9 (MATLAB 1..10).
#
# The interpolation loop uses x1 in range(19) and x2 in range(1,6)
# corresponding to MATLAB x1=1:19, x2=1:5.  The formula is identical.
#
# MATLAB's  handles.Variables  is stored as  handles['Variables']  (a dict).
# MATLAB structs become Python dicts, consistent with NumberCrunching_100000.py
# and Evaluation.py.
# ---------------------------------------------------------------------------


# ===================================================================
#  COLONOSCOPY LIKELYHOOD COMPUTATION (core logic)
# ===================================================================

def compute_colonoscopy_likelyhood(colonoscopy_rate):
    """
    Compute the interpolated ColonoscopyLikelyhood array from ColonoscopyRate.

    This replicates the MATLAB logic:
        counter = 1;
        for x1=1:19
            for x2=1:5
                ColonoscopyLikelyhood(counter) = (ColonoscopyRate(x1) * (5-x2) + ...
                    ColonoscopyRate(x1+1) * (x2-1)) / 4;
                counter = counter + 1;
            end
        end
        ColonoscopyLikelyhood(counter : 150) = ColonoscopyLikelyhood(end);

    Parameters
    ----------
    colonoscopy_rate : array-like
        Array of 20 colonoscopy rate values (0-based indices 0..19).

    Returns
    -------
    numpy.ndarray
        Array of 150 interpolated colonoscopy likelyhood values (0-based
        indices 0..149).
    """
    colonoscopy_rate = np.asarray(colonoscopy_rate, dtype=float)
    colonoscopy_likelyhood = np.zeros(150)
    counter = 0  # 0-based (MATLAB counter starts at 1)
    for x1 in range(19):  # MATLAB x1=1:19 -> Python 0..18
        for x2 in range(1, 6):  # MATLAB x2=1:5 -> Python 1..5
            # MATLAB: (ColonoscopyRate(x1) * (5-x2) + ColonoscopyRate(x1+1) * (x2-1)) / 4
            # Python: colonoscopy_rate[x1] * (5 - x2) + colonoscopy_rate[x1 + 1] * (x2 - 1)) / 4
            colonoscopy_likelyhood[counter] = (
                colonoscopy_rate[x1] * (5 - x2) +
                colonoscopy_rate[x1 + 1] * (x2 - 1)
            ) / 4.0
            counter += 1
    # counter is now 95 (19*5). MATLAB: ColonoscopyLikelyhood(counter:150) = end value
    # In MATLAB, 'end' refers to the last assigned element (counter-1 in MATLAB = index 94 in Python
    # which is counter-1=94 in 0-based). Fill indices 95..149 with the value at index 94.
    colonoscopy_likelyhood[counter:150] = colonoscopy_likelyhood[counter - 1]
    return colonoscopy_likelyhood


# ===================================================================
#  SAFE NUMERIC PARSING  (replaces MATLAB str2num)
# ===================================================================

def _safe_str2num(text):
    """
    Parse a string to a float.  Returns (value, True) on success,
    or (None, False) on failure.  Mirrors MATLAB's [num, succ] = str2num(tmp).
    """
    try:
        value = float(text)
        return value, True
    except (ValueError, TypeError):
        return None, False


# ===================================================================
#  COLONOSCOPY SETTINGS GUI
# ===================================================================

class ColonoscopySettingsGUI:
    """
    Python/tkinter equivalent of the MATLAB GUIDE Colonoscopy_settings dialog.

    The GUI allows the user to view and edit:
      - ColonoscopyRate (20 values, displayed in text entries)
      - Colo_Detection (10 values)
      - RectoSigmo_Detection (10 values)
      - Risk parameters: perforation, serosa burn, bleeding, severe bleeding
      - Death parameters: perforation, severe bleeding
      - Rectosigmo perforation risk

    Changes are applied to the Variables dict and the interpolated
    ColonoscopyLikelyhood array is recomputed live.  The user may accept or
    reject the changes via the Done button.

    Parameters
    ----------
    variables : dict
        The simulation Variables dictionary (handles['Variables']).
        This dictionary is modified in-place if the user accepts changes.
    parent : tk.Tk or tk.Toplevel, optional
        Parent window.  If None, a new Tk root is created.
    """

    # Internal names for the 20 ColonoscopyRate entry widgets
    COLONOSCOPY_RATE_HANDLES = [
        'rate_1', 'rate_6', 'rate_11', 'rate_16',
        'rate_21', 'rate_26', 'rate_31', 'rate_36', 'rate_41',
        'rate_46', 'rate_51', 'rate_56', 'rate_61',
        'rate_66', 'rate_71', 'rate_76', 'rate_81',
        'rate_86', 'rate_91', 'rate_96',
    ]

    # Internal names for the 10 Colo_Detection entry widgets
    COLO_DETECTION_HANDLES = [
        'colo_detection_P1', 'colo_detection_P2', 'colo_detection_P3',
        'colo_detection_P4', 'colo_detection_P5', 'colo_detection_cis',
        'colo_detection_I', 'colo_detection_II', 'colo_detection_III',
        'colo_detection_IV',
    ]

    # Internal names for the 10 RectoSigmo_Detection entry widgets
    RECTOSIGMO_DETECTION_HANDLES = [
        'rectosigmo_detection_P1', 'rectosigmo_detection_P2',
        'rectosigmo_detection_P3', 'rectosigmo_detection_P4',
        'rectosigmo_detection_P5', 'rectosigmo_detection_P6',
        'rectosigmo_detection_Ca1', 'rectosigmo_detection_Ca2',
        'rectosigmo_detection_Ca3', 'rectosigmo_detection_Ca4',
    ]

    def __init__(self, variables, parent=None):
        # Store working copy and backup for cancel/revert
        self.variables = variables
        self.old_variables = copy.deepcopy(variables)

        # --- Build the window ------------------------------------------------
        if parent is None:
            self.root = tk.Tk()
        else:
            self.root = tk.Toplevel(parent)

        self.root.title('Colonoscopy')
        self.root.configure(bg='#9999ff')  # [0.6 0.6 1] in MATLAB

        # Dictionary to hold all widget references (mirrors MATLAB handles.xxx)
        self.widgets = {}

        self._build_gui()
        self._make_images_current()

    # -----------------------------------------------------------------
    #  GUI construction
    # -----------------------------------------------------------------
    def _build_gui(self):
        """Build all widgets."""
        root = self.root

        # --- Colonoscopy Rate section ----------------------------------------
        rate_frame = tk.LabelFrame(root, text='Colonoscopy Rate (ages 1-96 in steps of 5)',
                                   bg='#9999ff')
        rate_frame.pack(padx=5, pady=5, fill='x')

        age_labels = [
            '1', '6', '11', '16', '21', '26', '31', '36', '41',
            '46', '51', '56', '61', '66', '71', '76', '81',
            '86', '91', '96',
        ]
        for idx, name in enumerate(self.COLONOSCOPY_RATE_HANDLES):
            row = idx // 5
            col = idx % 5
            sub = tk.Frame(rate_frame, bg='#9999ff')
            sub.grid(row=row, column=col, padx=2, pady=2)
            tk.Label(sub, text=f'Age {age_labels[idx]}:', bg='#9999ff').pack(side='left')
            entry = tk.Entry(sub, width=8, state='disabled')
            entry.pack(side='left')
            self.widgets[name] = entry

        # --- Colo Detection section ------------------------------------------
        colo_frame = tk.LabelFrame(root, text='Colonoscopy Detection Rates',
                                   bg='#9999ff')
        colo_frame.pack(padx=5, pady=5, fill='x')

        colo_labels = ['P1', 'P2', 'P3', 'P4', 'P5', 'CIS', 'I', 'II', 'III', 'IV']
        for idx, name in enumerate(self.COLO_DETECTION_HANDLES):
            row = idx // 5
            col = idx % 5
            sub = tk.Frame(colo_frame, bg='#9999ff')
            sub.grid(row=row, column=col, padx=2, pady=2)
            tk.Label(sub, text=f'{colo_labels[idx]}:', bg='#9999ff').pack(side='left')
            entry = tk.Entry(sub, width=8)
            entry.pack(side='left')
            # Bind callback: on Return key, update the corresponding detection value
            entry.bind('<Return>', lambda e, c=idx: self._colo_detection_callback(c))
            self.widgets[name] = entry

        # --- RectoSigmo Detection section ------------------------------------
        recto_frame = tk.LabelFrame(root, text='Rectosigmoidoscopy Detection Rates',
                                    bg='#9999ff')
        recto_frame.pack(padx=5, pady=5, fill='x')

        recto_labels = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'Ca1', 'Ca2', 'Ca3', 'Ca4']
        for idx, name in enumerate(self.RECTOSIGMO_DETECTION_HANDLES):
            row = idx // 5
            col = idx % 5
            sub = tk.Frame(recto_frame, bg='#9999ff')
            sub.grid(row=row, column=col, padx=2, pady=2)
            tk.Label(sub, text=f'{recto_labels[idx]}:', bg='#9999ff').pack(side='left')
            entry = tk.Entry(sub, width=8)
            entry.pack(side='left')
            entry.bind('<Return>', lambda e, c=idx: self._rectosigmo_detection_callback(c))
            self.widgets[name] = entry

        # --- Risks section ---------------------------------------------------
        risk_frame = tk.LabelFrame(root, text='Risks and Complications',
                                   bg='#9999ff')
        risk_frame.pack(padx=5, pady=5, fill='x')

        risk_entries = [
            ('Perforation risk:', 'risc_perforation'),
            ('Serosa burn risk:', 'risc_serosa_burn'),
            ('Bleeding risk:', 'risc_bleeding'),
            ('Severe bleeding risk:', 'risc_severe_bleeding'),
            ('Death from perforation:', 'death_perforation'),
            ('Death from severe bleeding:', 'death_severe_bleeding'),
            ('Rectosigmo perforation risk:', 'risc_perforation_rectosigmo'),
        ]
        for row_idx, (label_text, widget_name) in enumerate(risk_entries):
            tk.Label(risk_frame, text=label_text, bg='#9999ff').grid(
                row=row_idx, column=0, sticky='e', padx=2, pady=2)
            entry = tk.Entry(risk_frame, width=12)
            entry.grid(row=row_idx, column=1, padx=2, pady=2)
            entry.bind('<Return>', lambda e, n=widget_name: self._risk_callback(n))
            self.widgets[widget_name] = entry

        # --- Done button -----------------------------------------------------
        btn_frame = tk.Frame(root, bg='#9999ff')
        btn_frame.pack(padx=5, pady=10)
        tk.Button(btn_frame, text='Done', command=self._done_callback).pack()

    # -----------------------------------------------------------------
    #  Make Images Current  (update all GUI fields from Variables)
    # -----------------------------------------------------------------
    def _make_images_current(self):
        """
        Refresh all GUI widgets from self.variables and recompute
        ColonoscopyLikelyhood.

        This is the Python equivalent of the MATLAB function
        MakeImagesCurrent(hObject, handles).
        """
        variables = self.variables

        # --- Colonoscopy Rate entries (disabled, display-only) ---------------
        colonoscopy_rate = np.asarray(variables['ColonoscopyRate'], dtype=float)
        for f in range(len(self.COLONOSCOPY_RATE_HANDLES)):
            entry = self.widgets[self.COLONOSCOPY_RATE_HANDLES[f]]
            entry.config(state='normal')
            entry.delete(0, tk.END)
            entry.insert(0, str(colonoscopy_rate[f]))
            entry.config(state='disabled')

        # --- Compute ColonoscopyLikelyhood -----------------------------------
        # MATLAB:
        #   counter = 1;
        #   for x1=1:19
        #       for x2=1:5
        #           ColonoscopyLikelyhood(counter) = (ColonoscopyRate(x1)*(5-x2) +
        #               ColonoscopyRate(x1+1)*(x2-1)) / 4;
        #           counter = counter + 1;
        #       end
        #   end
        #   ColonoscopyLikelyhood(counter:150) = ColonoscopyLikelyhood(end);
        variables['ColonoscopyLikelyhood'] = list(
            compute_colonoscopy_likelyhood(colonoscopy_rate)
        )

        # --- Colo Detection entries ------------------------------------------
        colo_detection = variables['Colo_Detection']
        for f in range(len(self.COLO_DETECTION_HANDLES)):
            entry = self.widgets[self.COLO_DETECTION_HANDLES[f]]
            entry.delete(0, tk.END)
            entry.insert(0, str(colo_detection[f]))

        # --- RectoSigmo Detection entries ------------------------------------
        rectosigmo_detection = variables['RectoSigmo_Detection']
        for f in range(len(self.RECTOSIGMO_DETECTION_HANDLES)):
            entry = self.widgets[self.RECTOSIGMO_DETECTION_HANDLES[f]]
            entry.delete(0, tk.END)
            entry.insert(0, str(rectosigmo_detection[f]))

        # --- Risk entries ----------------------------------------------------
        self._set_entry('risc_perforation',
                        str(variables['Colonoscopy_RiscPerforation']))
        self._set_entry('risc_serosa_burn',
                        str(variables['Colonoscopy_RiscSerosaBurn']))
        self._set_entry('risc_bleeding',
                        str(variables['Colonoscopy_RiscBleeding']))
        self._set_entry('risc_severe_bleeding',
                        str(variables['Colonoscopy_RiscBleedingTransfusion']))
        self._set_entry('death_perforation',
                        str(variables['DeathPerforation']))
        self._set_entry('death_severe_bleeding',
                        str(variables['DeathBleedingTransfusion']))
        self._set_entry('risc_perforation_rectosigmo',
                        str(variables['Rectosigmo_Perforation']))

    def _set_entry(self, widget_name, text):
        """Helper: set the text of an Entry widget."""
        entry = self.widgets[widget_name]
        entry.delete(0, tk.END)
        entry.insert(0, text)

    def _get_entry(self, widget_name):
        """Helper: get the text content of an Entry widget."""
        return self.widgets[widget_name].get()

    # -----------------------------------------------------------------
    #  Rate callbacks  (one per ColonoscopyRate entry)
    # -----------------------------------------------------------------
    # In MATLAB, each rate_XX_Callback reads the entry, validates, and
    # updates Variables.ColonoscopyRate(c).  Here we generalize with a
    # single method.

    def _rate_callback(self, c):
        """
        Callback for ColonoscopyRate entry at 0-based index c.

        MATLAB equivalent (example for c=0, which is MATLAB c=1):
            c=1; tmp=get(handles.(handles.ColonoscopyRateHandles{c}), 'string');
            [num, succ] = str2num(tmp);
            if succ, if num >=0, handles.Variables.ColonoscopyRate(c)=num; end, end
            handles = MakeImagesCurrent(hObject, handles);
        """
        tmp = self._get_entry(self.COLONOSCOPY_RATE_HANDLES[c])
        num, succ = _safe_str2num(tmp)
        if succ:
            if num >= 0:
                self.variables['ColonoscopyRate'][c] = num
        self._make_images_current()

    # -----------------------------------------------------------------
    #  Colo Detection callbacks
    # -----------------------------------------------------------------

    def _colo_detection_callback(self, c):
        """
        Callback for Colo_Detection entry at 0-based index c.

        MATLAB equivalent (example for c=0, which is MATLAB c=1):
            c=1; tmp=get(handles.(handles.ColoDetectionHandles{c}), 'string');
            [num, succ] = str2num(tmp);
            if succ, if num >=0, handles.Variables.Colo_Detection(c)=num; end, end
            handles = MakeImagesCurrent(hObject, handles);
        """
        tmp = self._get_entry(self.COLO_DETECTION_HANDLES[c])
        num, succ = _safe_str2num(tmp)
        if succ:
            if num >= 0:
                self.variables['Colo_Detection'][c] = num
        self._make_images_current()

    # -----------------------------------------------------------------
    #  RectoSigmo Detection callbacks
    # -----------------------------------------------------------------

    def _rectosigmo_detection_callback(self, c):
        """
        Callback for RectoSigmo_Detection entry at 0-based index c.

        MATLAB equivalent (example for c=0, which is MATLAB c=1):
            c=1; tmp=get(handles.(handles.RectosigmoDetectionHandles{c}), 'string');
            [num, succ] = str2num(tmp);
            if succ, if num >=0, handles.Variables.RectoSigmo_Detection(c)=num; end, end
            handles = MakeImagesCurrent(hObject, handles);
        """
        tmp = self._get_entry(self.RECTOSIGMO_DETECTION_HANDLES[c])
        num, succ = _safe_str2num(tmp)
        if succ:
            if num >= 0:
                self.variables['RectoSigmo_Detection'][c] = num
        self._make_images_current()

    # -----------------------------------------------------------------
    #  Risk / complication callbacks
    # -----------------------------------------------------------------

    def _risk_callback(self, widget_name):
        """
        Callback for risk/complication entries.

        Each callback reads the entry value, validates it, and updates the
        corresponding Variables key.

        MATLAB equivalents:
            risc_perforation_Callback     -> Variables.Colonoscopy_Risc_Normal
            risc_serosa_burn_Callback     -> Variables.Colonoscopy_Death_Normal
            risc_bleeding_Callback        -> Variables.Colonoscopy_Risc_Polyp
            risc_severe_bleeding_Callback -> Variables.Colonoscopy_Death_Polyp
            death_perforation_Callback    -> Variables.DeathPerforation
            death_severe_bleeding_Callback-> Variables.DeathBleedingTransfusion
            risc_perforation_rectosigmo_Callback -> Variables.Rectosigmo_Perforation

        Note: The MATLAB callbacks update different variable names than the
        display values shown.  The display shows Colonoscopy_RiscPerforation
        etc., but the callbacks write to Colonoscopy_Risc_Normal, etc.
        This matches the original MATLAB code exactly.
        """
        # Map widget name to the Variables key that the MATLAB callback updates
        callback_var_map = {
            'risc_perforation': 'Colonoscopy_Risc_Normal',
            'risc_serosa_burn': 'Colonoscopy_Death_Normal',
            'risc_bleeding': 'Colonoscopy_Risc_Polyp',
            'risc_severe_bleeding': 'Colonoscopy_Death_Polyp',
            'death_perforation': 'DeathPerforation',
            'death_severe_bleeding': 'DeathBleedingTransfusion',
            'risc_perforation_rectosigmo': 'Rectosigmo_Perforation',
        }

        var_key = callback_var_map[widget_name]
        tmp = self._get_entry(widget_name)
        num, succ = _safe_str2num(tmp)
        if succ:
            if num >= 0:
                self.variables[var_key] = num
        self._make_images_current()

    # -----------------------------------------------------------------
    #  Done callback
    # -----------------------------------------------------------------

    def _done_callback(self):
        """
        Done button callback.

        MATLAB equivalent:
            Answer = questdlg('Do you want to keep the settings?', ...
                              'Return?', 'Yes', 'No', 'Cancel', 'Yes');
            if isequal(Answer, 'Cancel'), return, end
            if isequal(Answer, 'No')
                handles.Variables = handles.OldVariables;
            end
            set(0, 'userdata', handles.Variables);
            uiresume(handles.figure1);
            if ishandle(handles.figure1), delete(handles.figure1); end
        """
        answer = messagebox.askyesnocancel(
            'Return?', 'Do you want to keep the settings?')
        if answer is None:
            # Cancel
            return
        elif answer is False:
            # No -- revert to old variables
            # Restore all keys from the backup
            self.variables.clear()
            self.variables.update(copy.deepcopy(self.old_variables))
        # else: Yes -- keep current variables (already in self.variables)

        self.root.destroy()

    # -----------------------------------------------------------------
    #  Run  (blocking main loop)
    # -----------------------------------------------------------------
    def run(self):
        """Start the tkinter main loop (blocking)."""
        self.root.mainloop()


# ===================================================================
#  MODULE-LEVEL CONVENIENCE FUNCTION
# ===================================================================

def colonoscopy_settings(variables):
    """
    Open the Colonoscopy settings dialog and return the (possibly modified)
    Variables dictionary.

    This is the main entry point, equivalent to calling the MATLAB function
    ``Colonoscopy_settings`` which opens a GUIDE dialog, blocks until the
    user clicks Done, and returns the updated userdata.

    Parameters
    ----------
    variables : dict
        The simulation Variables dictionary (handles['Variables']).
        Modified in-place if the user accepts changes.

    Returns
    -------
    dict
        The same ``variables`` dict reference (possibly reverted if the user
        chose 'No').
    """
    gui = ColonoscopySettingsGUI(variables)
    gui.run()
    return variables


# ===================================================================
#  STANDALONE EXECUTION
# ===================================================================

if __name__ == '__main__':
    # Example: create a minimal Variables dict for testing
    test_variables = {
        'ColonoscopyRate': [0.0] * 20,
        'ColonoscopyLikelyhood': [0.0] * 150,
        'Colo_Detection': [0.65, 0.75, 0.81, 0.87, 0.95,
                           0.95, 0.95, 0.95, 0.95, 1.0],
        'RectoSigmo_Detection': [0.4, 0.5, 0.56, 0.62, 0.825,
                                 0.825, 0.95, 0.95, 0.95, 1.0],
        'Colonoscopy_RiscPerforation': 0.0007,
        'Colonoscopy_RiscSerosaBurn': 0.0003,
        'Colonoscopy_RiscBleeding': 0.0011,
        'Colonoscopy_RiscBleedingTransfusion': 0.0004,
        'DeathPerforation': 0.052,
        'DeathBleedingTransfusion': 0.0052,
        'Rectosigmo_Perforation': 2e-05,
        # Additional keys that the risk callbacks write to
        'Colonoscopy_Risc_Normal': 0.0007,
        'Colonoscopy_Death_Normal': 0.0003,
        'Colonoscopy_Risc_Polyp': 0.0011,
        'Colonoscopy_Death_Polyp': 0.0004,
    }
    result = colonoscopy_settings(test_variables)
    print('Final ColonoscopyLikelyhood (first 10):', result['ColonoscopyLikelyhood'][:10])
    print('Final Colo_Detection:', result['Colo_Detection'])
