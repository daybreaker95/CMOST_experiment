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
Mortality_Settings.py -- Python translation of Mortality_Settings.m

Original MATLAB purpose:
    GUIDE-based GUI for configuring mortality correction parameters.
    Allows editing of 20 MortalityCorrection values (at age intervals of 5).
    Computes and displays the interpolated MortalityCorrectionGraph (150 values)
    with a matplotlib plot.

Parameters (via constructor):
    variables : dict
        The simulation Variables dictionary containing 'MortalityCorrection'
        (list of 20 values) and 'MortalityCorrectionGraph' (list of 150 values).

Returns (via module-level function):
    dict : The same variables dict reference (possibly reverted).
"""

import numpy as np
import copy
import tkinter as tk
from tkinter import messagebox

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

# ---------------------------------------------------------------------------
# INDEX CONVENTION NOTES:
#
# MortalityCorrection has 20 elements: index 0..19 (MATLAB 1..20).
# MortalityCorrectionGraph has 150 elements: index 0..149 (MATLAB 1..150).
# The interpolation logic is identical to Colonoscopy_settings.py's
# compute_colonoscopy_likelyhood pattern.
# ---------------------------------------------------------------------------


def compute_mortality_correction_graph(mortality_correction):
    """
    Compute the interpolated MortalityCorrectionGraph from MortalityCorrection.

    MATLAB logic:
        counter = 1;
        for x1=1:19
            for x2=1:5
                MortalityCorrectionGraph(counter) = (MortalityCorrection(x1)*(5-x2) +
                    MortalityCorrection(x1+1)*(x2-1))/4;
                counter = counter + 1;
            end
        end
        MortalityCorrectionGraph(counter:150) = MortalityCorrection(end);

    Parameters
    ----------
    mortality_correction : array-like
        20-element array of mortality correction values (0-based indices 0..19).

    Returns
    -------
    numpy.ndarray
        150-element interpolated array.
    """
    mc = np.asarray(mortality_correction, dtype=float)
    graph = np.zeros(150)
    counter = 0
    for x1 in range(19):  # MATLAB x1=1:19
        for x2 in range(1, 6):  # MATLAB x2=1:5
            graph[counter] = (mc[x1] * (5 - x2) + mc[x1 + 1] * (x2 - 1)) / 4.0
            counter += 1
    # Fill remaining with last MortalityCorrection value
    graph[counter:150] = mc[-1]
    return graph


def _safe_str2num(text):
    """Parse a string to float. Returns (value, True) on success."""
    try:
        return float(text), True
    except (ValueError, TypeError):
        return None, False


# Widget handle names matching MATLAB handles.MortHandles
MORT_HANDLES = [
    'Mort_01', 'Mort_06', 'Mort_11', 'Mort_16',
    'Mort_21', 'Mort_26', 'Mort_31', 'Mort_36', 'Mort_41',
    'Mort_46', 'Mort_51', 'Mort_56', 'Mort_61',
    'Mort_66', 'Mort_71', 'Mort_76', 'Mort_81',
    'Mort_86', 'Mort_91', 'Mort_96',
]


class MortalitySettingsGUI:
    """
    Python/tkinter equivalent of the MATLAB GUIDE Mortality_Settings dialog.

    Displays 20 editable mortality correction values and a plot of the
    interpolated 150-element MortalityCorrectionGraph.

    Parameters
    ----------
    variables : dict
        Simulation Variables dict. Must contain 'MortalityCorrection' (20 values).
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

        self.root.title('Mortality')
        self.root.configure(bg='#9999ff')
        self.widgets = {}

        self._build_gui()
        self._make_images_current()

    def _build_gui(self):
        root = self.root

        # Mortality correction entries
        entry_frame = tk.LabelFrame(root, text='Mortality Correction (ages 1-96 in steps of 5)',
                                    bg='#9999ff')
        entry_frame.pack(padx=5, pady=5, fill='x')

        age_labels = ['1', '6', '11', '16', '21', '26', '31', '36', '41',
                       '46', '51', '56', '61', '66', '71', '76', '81',
                       '86', '91', '96']

        for idx, name in enumerate(MORT_HANDLES):
            row = idx // 5
            col = idx % 5
            sub = tk.Frame(entry_frame, bg='#9999ff')
            sub.grid(row=row, column=col, padx=2, pady=2)
            tk.Label(sub, text=f'Age {age_labels[idx]}:', bg='#9999ff').pack(side='left')
            entry = tk.Entry(sub, width=8)
            entry.pack(side='left')
            entry.bind('<Return>', lambda e, c=idx: self._mort_callback(c))
            self.widgets[name] = entry

        # Plot area
        # MATLAB: axes(handles.Axis_1); plot(1:100, MortalityCorrectionGraph(1:100))
        self.fig = Figure(figsize=(6, 2.5), dpi=80)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().pack(padx=5, pady=5, fill='x')

        # Done button
        btn_frame = tk.Frame(root, bg='#9999ff')
        btn_frame.pack(padx=5, pady=10)
        tk.Button(btn_frame, text='Done', command=self._done_callback).pack()

    def _make_images_current(self):
        """
        Refresh entries, recompute MortalityCorrectionGraph, and update plot.

        MATLAB equivalent: MakeImagesCurrent(hObject, handles)
        """
        mc = self.variables.get('MortalityCorrection', [1.0] * 20)

        # Recompute the interpolated graph
        graph = compute_mortality_correction_graph(mc)
        self.variables['MortalityCorrectionGraph'] = list(graph)

        # Update entries
        for f in range(len(MORT_HANDLES)):
            entry = self.widgets[MORT_HANDLES[f]]
            entry.delete(0, tk.END)
            if f < len(mc):
                entry.insert(0, str(mc[f]))
            else:
                entry.insert(0, '1.0')

        # Update plot
        self.ax.clear()
        self.ax.plot(range(1, 101), graph[:100])
        # Plot marker points at every 5th year
        for f in range(min(20, len(mc))):
            self.ax.plot(f * 5 + 1, mc[f], 'gs', markersize=3)
        self.ax.set_title('Mortality Correction Graph')
        self.canvas.draw()

    def _mort_callback(self, c):
        """
        Callback for MortalityCorrection entry at 0-based index c.

        MATLAB equivalent:
            c=<val>; tmp=get(handles.(...), 'string'); [num,succ]=str2num(tmp);
            if succ, if num>=0, handles.Variables.MortalityCorrection(c)=num; end, end
        """
        tmp = self.widgets[MORT_HANDLES[c]].get()
        num, succ = _safe_str2num(tmp)
        if succ and num >= 0:
            mc = self.variables.get('MortalityCorrection', [1.0] * 20)
            while len(mc) <= c:
                mc.append(1.0)
            mc[c] = num
            self.variables['MortalityCorrection'] = mc
        self._make_images_current()

    def _done_callback(self):
        """
        MATLAB equivalent: done_Callback
        """
        answer = messagebox.askyesnocancel('Return?', 'Do you want to keep the settings?')
        if answer is None:
            return
        elif answer is False:
            self.variables.clear()
            self.variables.update(copy.deepcopy(self.old_variables))
        self.root.destroy()

    def run(self):
        self.root.mainloop()


def mortality_settings(variables):
    """
    Open the Mortality Settings dialog.

    Parameters
    ----------
    variables : dict
        Simulation Variables dict.

    Returns
    -------
    dict
        The same variables dict reference.
    """
    gui = MortalitySettingsGUI(variables)
    gui.run()
    return variables


if __name__ == '__main__':
    test_vars = {
        'MortalityCorrection': [1.0] * 20,
        'MortalityCorrectionGraph': [1.0] * 150,
    }
    result = mortality_settings(test_vars)
    print('MortalityCorrection:', result['MortalityCorrection'])
