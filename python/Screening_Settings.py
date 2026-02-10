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

"""
Screening_Settings.py -- Python translation of Screening_Settings.m

Original MATLAB purpose:
    GUIDE-based GUI for configuring screening test parameters in CMOST.
    Allows editing of colonoscopy, rectosigmoidoscopy, FOBT, I-FOBT,
    Sept9 (high sensitivity / high specificity), and other screening
    test parameters including participation rates, adherence, follow-up,
    age windows, intervals, and test sensitivities.

Parameters (via constructor):
    variables : dict
        The simulation Variables dictionary containing a 'Screening' sub-dict.
        Modified in-place if the user accepts changes.

Returns (via module-level function):
    dict : The same variables dict reference (possibly reverted if user chose 'No').
"""

import copy
import tkinter as tk
from tkinter import messagebox


# ---------------------------------------------------------------------------
# INDEX CONVENTION NOTES:
#
# In the original MATLAB code, arrays are 1-based.  In this Python
# translation, all array indices are 0-based.
#
# Each screening test has 7 or 8 parameter values stored as a list.
# Colonoscopy has 7 values (no follow-up parameter).
# Other tests have 8 values each.
# Sensitivity arrays have 10 values each.
#
# MATLAB's handles.Variables.Screening.Colonoscopy(c) becomes
# self.variables['Screening']['Colonoscopy'][c-1] in Python.
# ---------------------------------------------------------------------------


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


# Screening test configuration: test name, variable key, parameter labels, has sensitivity
SCREENING_TESTS = [
    {
        'name': 'Colonoscopy',
        'key': 'Colonoscopy',
        'labels': ['% Pop', 'Adherence', 'Y Start', 'Y End', 'Interval', 'Y After Colo', 'Specificity'],
        'has_sensitivity': False,
    },
    {
        'name': 'Rectosigmoidoscopy',
        'key': 'Rectosigmoidoscopy',
        'labels': ['% Pop', 'Adherence', 'Follow Up', 'Y Start', 'Y End', 'Interval', 'Y After Colo', 'Specificity'],
        'has_sensitivity': False,
    },
    {
        'name': 'FOBT',
        'key': 'FOBT',
        'labels': ['% Pop', 'Adherence', 'Follow Up', 'Y Start', 'Y End', 'Interval', 'Y After Colo', 'Specificity'],
        'has_sensitivity': True,
        'sens_key': 'FOBT_Sens',
    },
    {
        'name': 'I-FOBT',
        'key': 'I_FOBT',
        'labels': ['% Pop', 'Adherence', 'Follow Up', 'Y Start', 'Y End', 'Interval', 'Y After Colo', 'Specificity'],
        'has_sensitivity': True,
        'sens_key': 'I_FOBT_Sens',
    },
    {
        'name': 'Sept9 Hi-Sens',
        'key': 'Sept9_HiSens',
        'labels': ['% Pop', 'Adherence', 'Follow Up', 'Y Start', 'Y End', 'Interval', 'Y After Colo', 'Specificity'],
        'has_sensitivity': True,
        'sens_key': 'Sept9_HiSens_Sens',
    },
    {
        'name': 'Sept9 Hi-Spec',
        'key': 'Sept9_HiSpec',
        'labels': ['% Pop', 'Adherence', 'Follow Up', 'Y Start', 'Y End', 'Interval', 'Y After Colo', 'Specificity'],
        'has_sensitivity': True,
        'sens_key': 'Sept9_HiSpec_Sens',
    },
    {
        'name': 'Other',
        'key': 'other',
        'labels': ['% Pop', 'Adherence', 'Follow Up', 'Y Start', 'Y End', 'Interval', 'Y After Colo', 'Specificity'],
        'has_sensitivity': True,
        'sens_key': 'other_Sens',
    },
]

SENS_LABELS = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'Ca1', 'Ca2', 'Ca3', 'Ca4']


class ScreeningSettingsGUI:
    """
    Python/tkinter equivalent of the MATLAB GUIDE Screening_Settings dialog.

    The GUI allows the user to view and edit all screening test parameters
    including participation, adherence, follow-up, age start/end, interval,
    years after colonoscopy, specificity, and test sensitivities for
    FOBT, I-FOBT, Sept9, and other tests.

    Parameters
    ----------
    variables : dict
        The simulation Variables dictionary containing 'Screening' sub-dict.
        Modified in-place if user accepts changes.
    parent : tk.Tk or tk.Toplevel, optional
        Parent window.
    """

    def __init__(self, variables, parent=None):
        self.variables = variables
        self.old_variables = copy.deepcopy(variables)

        if parent is None:
            self.root = tk.Tk()
        else:
            self.root = tk.Toplevel(parent)

        self.root.title('Screening')
        self.root.configure(bg='#9999ff')
        self.root.minsize(800, 600)

        # Widget storage
        self.widgets = {}  # key -> Entry widget

        self._build_gui()
        self._make_images_current()

    def _build_gui(self):
        """Build all widgets for screening settings."""
        root = self.root

        # Create a canvas with scrollbar for the many fields
        canvas = tk.Canvas(root, bg='#9999ff')
        scrollbar = tk.Scrollbar(root, orient='vertical', command=canvas.yview)
        scrollable_frame = tk.Frame(canvas, bg='#9999ff')
        scrollable_frame.bind('<Configure>', lambda e: canvas.configure(scrollregion=canvas.bbox('all')))
        canvas.create_window((0, 0), window=scrollable_frame, anchor='nw')
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')

        # Build each screening test section
        for test_info in SCREENING_TESTS:
            test_name = test_info['name']
            test_key = test_info['key']
            labels = test_info['labels']

            # Parameter frame
            frame = tk.LabelFrame(scrollable_frame, text=test_name, bg='#9999ff')
            frame.pack(padx=5, pady=3, fill='x')

            row_frame = tk.Frame(frame, bg='#9999ff')
            row_frame.pack(fill='x', padx=3, pady=2)

            for idx, label in enumerate(labels):
                sub = tk.Frame(row_frame, bg='#9999ff')
                sub.pack(side='left', padx=2)
                tk.Label(sub, text=label, bg='#9999ff', font=('TkDefaultFont', 8)).pack()
                entry = tk.Entry(sub, width=8)
                entry.pack()
                widget_key = f'{test_key}_param_{idx}'
                entry.bind('<Return>', lambda e, k=test_key, i=idx: self._param_callback(k, i))
                self.widgets[widget_key] = entry

            # Sensitivity section if applicable
            if test_info.get('has_sensitivity'):
                sens_key = test_info['sens_key']
                sens_frame = tk.Frame(frame, bg='#9999ff')
                sens_frame.pack(fill='x', padx=3, pady=2)
                tk.Label(sens_frame, text='Sensitivity:', bg='#9999ff', font=('TkDefaultFont', 8, 'bold')).pack(side='left')

                for idx, label in enumerate(SENS_LABELS):
                    sub = tk.Frame(sens_frame, bg='#9999ff')
                    sub.pack(side='left', padx=2)
                    tk.Label(sub, text=label, bg='#9999ff', font=('TkDefaultFont', 7)).pack()
                    entry = tk.Entry(sub, width=6)
                    entry.pack()
                    widget_key = f'{sens_key}_sens_{idx}'
                    entry.bind('<Return>', lambda e, sk=sens_key, i=idx: self._sens_callback(sk, i))
                    self.widgets[widget_key] = entry

        # Percent no screening display
        no_screen_frame = tk.Frame(scrollable_frame, bg='#9999ff')
        no_screen_frame.pack(padx=5, pady=3, fill='x')
        tk.Label(no_screen_frame, text='% No Screening:', bg='#9999ff').pack(side='left')
        no_screen_entry = tk.Entry(no_screen_frame, width=10, state='disabled')
        no_screen_entry.pack(side='left', padx=5)
        self.widgets['percent_no_screening'] = no_screen_entry

        # Done button
        btn_frame = tk.Frame(scrollable_frame, bg='#9999ff')
        btn_frame.pack(padx=5, pady=10)
        tk.Button(btn_frame, text='Done', command=self._done_callback, width=15).pack()

    def _make_images_current(self):
        """
        Refresh all GUI widgets from self.variables['Screening'].

        MATLAB equivalent: MakeImagesCurrent(hObject, handles)
        """
        screening = self.variables.get('Screening', {})

        for test_info in SCREENING_TESTS:
            test_key = test_info['key']
            labels = test_info['labels']
            values = screening.get(test_key, [0] * len(labels))

            for idx in range(len(labels)):
                widget_key = f'{test_key}_param_{idx}'
                if widget_key in self.widgets:
                    entry = self.widgets[widget_key]
                    entry.delete(0, tk.END)
                    if idx < len(values):
                        entry.insert(0, str(values[idx]))
                    else:
                        entry.insert(0, '0')

            if test_info.get('has_sensitivity'):
                sens_key = test_info['sens_key']
                sens_values = screening.get(sens_key, [0] * 10)
                for idx in range(10):
                    widget_key = f'{sens_key}_sens_{idx}'
                    if widget_key in self.widgets:
                        entry = self.widgets[widget_key]
                        entry.delete(0, tk.END)
                        if idx < len(sens_values):
                            entry.insert(0, str(sens_values[idx]))
                        else:
                            entry.insert(0, '0')

        # Calculate % no screening
        # MATLAB: Summe = sum of first element of each test
        summe = 0
        for test_info in SCREENING_TESTS:
            vals = screening.get(test_info['key'], [0])
            if len(vals) > 0:
                summe += vals[0]

        no_screen_entry = self.widgets['percent_no_screening']
        no_screen_entry.config(state='normal')
        no_screen_entry.delete(0, tk.END)
        no_screen_entry.insert(0, str(1 - summe))
        no_screen_entry.config(state='disabled')

    def _param_callback(self, test_key, idx):
        """
        Callback for screening test parameter entry.

        MATLAB equivalent: e.g. Colo_PercPop_Callback, FOBT_adherence_Callback, etc.
        All follow the same pattern:
            c=<idx>; tmp=get(handles.(...), 'string'); [num, succ]=str2num(tmp);
            if succ, if num>=0, handles.Variables.Screening.<key>(c)=num; end, end
        """
        widget_key = f'{test_key}_param_{idx}'
        tmp = self.widgets[widget_key].get()
        num, succ = _safe_str2num(tmp)
        if succ and num >= 0:
            if 'Screening' not in self.variables:
                self.variables['Screening'] = {}
            screening = self.variables['Screening']
            if test_key not in screening:
                screening[test_key] = [0] * 8
            # Ensure the list is long enough
            while len(screening[test_key]) <= idx:
                screening[test_key].append(0)
            screening[test_key][idx] = num
        self._make_images_current()

    def _sens_callback(self, sens_key, idx):
        """
        Callback for sensitivity parameter entry.

        MATLAB equivalent: e.g. FOBT_Sens_P1_Callback, etc.
        All follow: c=<idx>; ... handles.Variables.Screening.<sens_key>(c)=num;
        """
        widget_key = f'{sens_key}_sens_{idx}'
        tmp = self.widgets[widget_key].get()
        num, succ = _safe_str2num(tmp)
        if succ and num >= 0:
            screening = self.variables.get('Screening', {})
            if sens_key not in screening:
                screening[sens_key] = [0] * 10
            while len(screening[sens_key]) <= idx:
                screening[sens_key].append(0)
            screening[sens_key][idx] = num
        self._make_images_current()

    def _done_callback(self):
        """
        Done button callback.

        MATLAB equivalent: done_Callback
            Answer = questdlg('Do you want to keep the settings?', ...);
            if 'Cancel' return; if 'No' revert; set(0,'userdata',...); uiresume
        """
        answer = messagebox.askyesnocancel('Return?', 'Do you want to keep the settings?')
        if answer is None:
            return
        elif answer is False:
            self.variables.clear()
            self.variables.update(copy.deepcopy(self.old_variables))
        self.root.destroy()

    def run(self):
        """Start the tkinter main loop (blocking)."""
        self.root.mainloop()


def screening_settings(variables):
    """
    Open the Screening Settings dialog and return the (possibly modified)
    Variables dictionary.

    Parameters
    ----------
    variables : dict
        The simulation Variables dictionary containing 'Screening' sub-dict.

    Returns
    -------
    dict
        The same variables dict reference (possibly reverted if user chose 'No').
    """
    gui = ScreeningSettingsGUI(variables)
    gui.run()
    return variables


if __name__ == '__main__':
    test_variables = {
        'Screening': {
            'Mode': 'off',
            'Colonoscopy': [0, 0.75, 50, 81, 10, 5, 0],
            'Rectosigmoidoscopy': [0, 0.75, 0.9, 50, 81, 5, 5, 0],
            'FOBT': [0, 0.5, 0.9, 50, 81, 1, 5, 0.98],
            'I_FOBT': [0, 0.5, 0.9, 50, 81, 1, 5, 0.95],
            'Sept9_HiSens': [0, 0.9, 0.9, 50, 81, 1, 5, 0.85],
            'Sept9_HiSpec': [0, 0.9, 0.9, 50, 81, 1, 5, 0.99],
            'other': [0, 0, 0, 50, 81, 1, 5, 0.925],
            'FOBT_Sens': [0.02, 0.02, 0.05, 0.05, 0.12, 0.12, 0.4, 0.4, 0.4, 0.4],
            'I_FOBT_Sens': [0.05, 0.05, 0.101, 0.101, 0.22, 0.22, 0.7, 0.7, 0.7, 0.7],
            'Sept9_HiSens_Sens': [0, 0, 0, 0, 0, 0, 0.89, 0.93, 0.99, 0.99],
            'Sept9_HiSpec_Sens': [0, 0, 0, 0, 0, 0, 0.67, 0.86, 0.87, 0.82],
            'other_Sens': [0.075, 0.075, 0.124, 0.124, 0.239, 0.239, 0.7, 0.7, 0.7, 0.7],
        }
    }
    result = screening_settings(test_variables)
    print('Screening settings:', result['Screening'])
