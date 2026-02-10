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
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

"""
Location_Settings.py -- Python translation of Location_Settings.m

Original MATLAB purpose:
    GUIDE-based GUI for configuring location-specific parameters in CMOST.
    Allows editing of 13-element arrays for:
      - New polyp location (Location_NewPolyp)
      - Early polyp progression (Location_EarlyProgression)
      - Advanced polyp progression (Location_AdvancedProgression)
      - Colonoscopy detection by location (Location_ColoDetection)
      - Rectosigmoidoscopy detection by location (Location_RectoSigmoDetection)
      - Colonoscopy reach by location (Location_ColoReach)
      - Rectosigmoidoscopy reach by location (Location_RectoSigmoReach)
      - Direct cancer location (Location_DirectCa)

Parameters (via constructor):
    variables : dict
        The simulation Variables dictionary with Location_* keys.

Returns (via module-level function):
    dict : The same variables dict reference (possibly reverted).
"""

import copy
import tkinter as tk
from tkinter import messagebox

# ---------------------------------------------------------------------------
# INDEX CONVENTION NOTES:
#
# Each location array has 13 elements: index 0..12 (MATLAB 1..13).
# The 13 locations represent different segments of the colon.
# ---------------------------------------------------------------------------


def _safe_str2num(text):
    """Parse a string to float. Returns (value, True) on success."""
    try:
        return float(text), True
    except (ValueError, TypeError):
        return None, False


# Configuration for the 8 location parameter groups
LOCATION_GROUPS = [
    {'label': 'New Polyp Location', 'var_key': 'Location_NewPolyp'},
    {'label': 'Early Progression', 'var_key': 'Location_EarlyProgression'},
    {'label': 'Advanced Progression', 'var_key': 'Location_AdvancedProgression'},
    {'label': 'Colonoscopy Detection', 'var_key': 'Location_ColoDetection'},
    {'label': 'Rectosigmo Detection', 'var_key': 'Location_RectoSigmoDetection'},
    {'label': 'Colonoscopy Reach', 'var_key': 'Location_ColoReach'},
    {'label': 'Rectosigmo Reach', 'var_key': 'Location_RectoSigmoReach'},
    {'label': 'Direct Cancer', 'var_key': 'Location_DirectCa'},
]


class LocationSettingsGUI:
    """
    Python/tkinter equivalent of the MATLAB GUIDE Location_Settings dialog.

    Displays 8 groups of 13 editable location-specific parameters.

    Parameters
    ----------
    variables : dict
        Simulation Variables dict with Location_* arrays.
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

        self.root.title('Location')
        self.root.configure(bg='#9999ff')
        self.root.minsize(750, 500)
        self.widgets = {}

        self._build_gui()
        self._make_images_current()

    def _build_gui(self):
        root = self.root

        # Scrollable area
        canvas = tk.Canvas(root, bg='#9999ff')
        scrollbar = tk.Scrollbar(root, orient='vertical', command=canvas.yview)
        scrollable_frame = tk.Frame(canvas, bg='#9999ff')
        scrollable_frame.bind('<Configure>', lambda e: canvas.configure(scrollregion=canvas.bbox('all')))
        canvas.create_window((0, 0), window=scrollable_frame, anchor='nw')
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')

        for group in LOCATION_GROUPS:
            var_key = group['var_key']
            frame = tk.LabelFrame(scrollable_frame, text=group['label'], bg='#9999ff')
            frame.pack(padx=5, pady=3, fill='x')

            row_frame = tk.Frame(frame, bg='#9999ff')
            row_frame.pack(fill='x', padx=3, pady=2)

            for idx in range(13):
                sub = tk.Frame(row_frame, bg='#9999ff')
                sub.pack(side='left', padx=1)
                tk.Label(sub, text=str(idx + 1), bg='#9999ff', font=('TkDefaultFont', 7)).pack()
                entry = tk.Entry(sub, width=6)
                entry.pack()
                widget_key = f'{var_key}_{idx}'
                entry.bind('<Return>', lambda e, vk=var_key, i=idx: self._value_callback(vk, i))
                self.widgets[widget_key] = entry

        # Done button
        btn_frame = tk.Frame(scrollable_frame, bg='#9999ff')
        btn_frame.pack(padx=5, pady=10)
        tk.Button(btn_frame, text='Done', command=self._done_callback, width=15).pack()

    def _make_images_current(self):
        """
        Refresh all entries from self.variables.

        MATLAB equivalent: MakeImagesCurrent
        """
        for group in LOCATION_GROUPS:
            var_key = group['var_key']
            values = self.variables.get(var_key, [0] * 13)
            for idx in range(13):
                widget_key = f'{var_key}_{idx}'
                if widget_key in self.widgets:
                    entry = self.widgets[widget_key]
                    entry.delete(0, tk.END)
                    if idx < len(values):
                        entry.insert(0, str(values[idx]))
                    else:
                        entry.insert(0, '0')

    def _value_callback(self, var_key, idx):
        """
        Generic callback for location parameter entry.

        MATLAB equivalent: new_X_Callback, early_X_Callback, etc.
        """
        widget_key = f'{var_key}_{idx}'
        tmp = self.widgets[widget_key].get()
        num, succ = _safe_str2num(tmp)
        if succ and num >= 0:
            values = self.variables.get(var_key, [0] * 13)
            if isinstance(values, list):
                while len(values) <= idx:
                    values.append(0)
                values[idx] = num
                self.variables[var_key] = values
            else:
                # numpy array or similar
                values[idx] = num
        self._make_images_current()

    def _done_callback(self):
        """MATLAB equivalent: done_Callback"""
        answer = messagebox.askyesnocancel('Return?', 'Do you want to keep the settings?')
        if answer is None:
            return
        elif answer is False:
            self.variables.clear()
            self.variables.update(copy.deepcopy(self.old_variables))
        self.root.destroy()

    def run(self):
        self.root.mainloop()


def location_settings(variables):
    """
    Open the Location Settings dialog.

    Parameters
    ----------
    variables : dict

    Returns
    -------
    dict
    """
    gui = LocationSettingsGUI(variables)
    gui.run()
    return variables


if __name__ == '__main__':
    test_vars = {
        'Location_NewPolyp': [1.0] * 13,
        'Location_EarlyProgression': [1.0] * 13,
        'Location_AdvancedProgression': [1.0] * 13,
        'Location_ColoDetection': [1.0] * 13,
        'Location_RectoSigmoDetection': [1.0] * 13,
        'Location_ColoReach': [1.0] * 13,
        'Location_RectoSigmoReach': [1.0] * 13,
        'Location_DirectCa': [1.0] * 13,
    }
    result = location_settings(test_vars)
    print('Location_NewPolyp:', result['Location_NewPolyp'])
