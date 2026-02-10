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

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

# ---------------------------------------------------------------------------
# INDEX CONVENTION NOTES:
#
# In the original MATLAB code, arrays are 1-based.  In this Python
# translation, all array indices are 0-based.
#
# Ind_Risk_Percentiles, Early_Risk_Percentiles, Adv_Risk_Percentiles each
# have 12 elements: index 0..11 (MATLAB 1..12).
#
# IndividualRisk, EarlyRisk, AdvRisk each have 500 elements: index 0..499
# (MATLAB 1..500).
#
# The percentile positions (Values) are [10,20,30,40,50,60,70,80,90,95,97,100]*5,
# which in MATLAB are 1-based indices into the 500-element arrays.
# In Python these become 0-based: we subtract 1 when indexing into arrays,
# but the AdjustRiskGraph computation is converted so that the resulting
# 500-element arrays are identical to the MATLAB output.
#
# MATLAB's  handles.Variables  is stored as  self.variables  (a dict).
# MATLAB structs become Python dicts, consistent with NumberCrunching_100000.py,
# Evaluation.py, Colonoscopy_settings.py, and Cost_Settings.py.
# ---------------------------------------------------------------------------


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
#  ADJUST RISK GRAPH  (core computation logic)
# ===================================================================

def adjust_risk_graph(variables):
    """
    Compute the interpolated risk arrays (IndividualRisk, EarlyRisk, AdvRisk)
    from the 12-element percentile arrays, and optionally normalize them.

    This replicates the MATLAB function AdjustRiskGraph(Variables) exactly.

    MATLAB logic (for IndividualRisk -- EarlyRisk and AdvRisk are identical
    in structure):

        Values = [10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 97, 100]*5;

        % individual polyp risk
        tmp(1:Values(1)) = Variables.Ind_Risk_Percentiles(1);
        for x1=1:length(Values)-1
            Start = Values(x1)+1;
            Ende  = Values(x1+1);
            for x2=Start:Ende
                tmp(x2) = (Variables.Ind_Risk_Percentiles(x1) * (Ende-x2) + ...
                    Variables.Ind_Risk_Percentiles(x1+1) * (x2-Start))/(Ende-Start);
            end
        end

        if isequal(Variables.RiskNormalize, 'on')
            tmp = tmp/mean(tmp);
            Variables.Ind_Risk_Percentiles = Variables.Ind_Risk_Percentiles/mean(tmp);
        end
        Variables.IndividualRisk = tmp;

    INDEX CONVERSION NOTES:
    -----------------------
    MATLAB Values(1) = 50.  tmp(1:50) means indices 1..50 in MATLAB,
    which is indices 0..49 in Python.

    MATLAB Start = Values(x1)+1, Ende = Values(x1+1).
    For x1=1 (MATLAB): Start=51, Ende=100; loop x2=51:100.
    In Python x1=0: Start=Values[0]+1=51, Ende=Values[1]=100.
    MATLAB tmp(x2) uses 1-based x2.  Python tmp[x2-1] uses 0-based.
    The interpolation formula references Ind_Risk_Percentiles(x1) and
    Ind_Risk_Percentiles(x1+1) in MATLAB (1-based), which become [x1]
    and [x1+1] in Python (0-based), since the MATLAB loop x1=1:11
    maps to Python x1=0:10.

    Parameters
    ----------
    variables : dict
        The simulation Variables dictionary.  Modified in-place:
        - variables['IndividualRisk'] is set to a 500-element list
        - variables['EarlyRisk'] is set to a 500-element list
        - variables['AdvRisk'] is set to a 500-element list
        - Percentile arrays may be modified if normalization is on.

    Returns
    -------
    dict
        The same ``variables`` dict reference.
    """
    # MATLAB: Values = [10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 97, 100]*5;
    Values = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 97, 100]) * 5
    # Values = [50, 100, 150, 200, 250, 300, 350, 400, 450, 475, 485, 500]

    # --- individual polyp risk ---
    ind_perc = np.asarray(variables['Ind_Risk_Percentiles'], dtype=float)
    tmp = np.zeros(500)
    # MATLAB: tmp(1:Values(1)) = Variables.Ind_Risk_Percentiles(1);
    # In MATLAB, Values(1) = 50, so tmp(1:50) = percentile(1).
    # In Python, Values[0] = 50, so tmp[0:50] = percentile[0].
    tmp[0:Values[0]] = ind_perc[0]

    # MATLAB: for x1=1:length(Values)-1  (x1 = 1..11)
    for x1 in range(len(Values) - 1):  # x1 = 0..10
        # MATLAB: Start = Values(x1)+1;  Ende = Values(x1+1);
        Start = Values[x1] + 1   # 1-based MATLAB index
        Ende = Values[x1 + 1]    # 1-based MATLAB index
        # MATLAB: for x2=Start:Ende
        for x2 in range(Start, Ende + 1):  # x2 in MATLAB 1-based range
            # MATLAB: tmp(x2) = (percentiles(x1) * (Ende-x2) +
            #                    percentiles(x1+1) * (x2-Start)) / (Ende-Start);
            # Python: tmp[x2-1] because x2 is 1-based
            tmp[x2 - 1] = (ind_perc[x1] * (Ende - x2) +
                           ind_perc[x1 + 1] * (x2 - Start)) / (Ende - Start)

    if variables.get('RiskNormalize', 'off') == 'on':
        mean_val = np.mean(tmp)
        if mean_val != 0:
            tmp = tmp / mean_val
            variables['Ind_Risk_Percentiles'] = list(
                np.asarray(variables['Ind_Risk_Percentiles'], dtype=float) / np.mean(tmp)
            )
    variables['IndividualRisk'] = list(tmp)

    # --- early progression risk ---
    early_perc = np.asarray(variables['Early_Risk_Percentiles'], dtype=float)
    tmp = np.zeros(500)
    tmp[0:Values[0]] = early_perc[0]

    for x1 in range(len(Values) - 1):
        Start = Values[x1] + 1
        Ende = Values[x1 + 1]
        for x2 in range(Start, Ende + 1):
            tmp[x2 - 1] = (early_perc[x1] * (Ende - x2) +
                           early_perc[x1 + 1] * (x2 - Start)) / (Ende - Start)

    if variables.get('RiskEarlyNormalize', 'off') == 'on':
        mean_val = np.mean(tmp)
        if mean_val != 0:
            tmp = tmp / mean_val
            variables['Early_Risk_Percentiles'] = list(
                np.asarray(variables['Early_Risk_Percentiles'], dtype=float) / np.mean(tmp)
            )
    variables['EarlyRisk'] = list(tmp)

    # --- advanced progression risk ---
    adv_perc = np.asarray(variables['Adv_Risk_Percentiles'], dtype=float)
    tmp = np.zeros(500)
    tmp[0:Values[0]] = adv_perc[0]

    for x1 in range(len(Values) - 1):
        Start = Values[x1] + 1
        Ende = Values[x1 + 1]
        for x2 in range(Start, Ende + 1):
            tmp[x2 - 1] = (adv_perc[x1] * (Ende - x2) +
                           adv_perc[x1 + 1] * (x2 - Start)) / (Ende - Start)

    if variables.get('RiskAdvNormalize', 'off') == 'on':
        mean_val = np.mean(tmp)
        if mean_val != 0:
            tmp = tmp / mean_val
            variables['Adv_Risk_Percentiles'] = list(
                np.asarray(variables['Adv_Risk_Percentiles'], dtype=float) / np.mean(tmp)
            )
    variables['AdvRisk'] = list(tmp)

    return variables


# ===================================================================
#  RISK SETTINGS GUI
# ===================================================================

class RiskSettingsGUI:
    """
    Python/tkinter equivalent of the MATLAB GUIDE Risk_Settings dialog.

    The GUI allows the user to view and edit risk percentile values for three
    risk categories:
      - Individual risk (Ind_Risk_Percentiles, 12 values)
      - Early progression risk (Early_Risk_Percentiles, 12 values)
      - Advanced progression risk (Adv_Risk_Percentiles, 12 values)

    Each category has:
      - A bar chart showing the interpolated 500-element risk array
      - Green markers at percentile positions
      - A red reference line at y=1
      - A log scale checkbox
      - A normalize checkbox

    Additionally there is a "correlate" checkbox that controls
    RiskCorrelation.

    Changes are applied to the Variables dict and the interpolated risk
    arrays are recomputed live.  The user may accept or reject the changes
    via the Done button.

    Parameters
    ----------
    variables : dict
        The simulation Variables dictionary (handles['Variables']).
        This dictionary is modified in-place if the user accepts changes.
    parent : tk.Tk or tk.Toplevel, optional
        Parent window.  If None, a new Tk root is created.
    """

    # Internal names for the 12 Individual Risk percentile entry widgets
    # (MATLAB: handles.IndividualRiskHandles)
    INDIVIDUAL_RISK_HANDLES = [
        'risk_10', 'risk_20', 'risk_30', 'risk_40',
        'risk_50', 'risk_60', 'risk_70', 'risk_80', 'risk_90',
        'risk_95', 'risk_97', 'risk_100',
    ]

    # Internal names for the 12 Early Risk percentile entry widgets
    # (MATLAB: handles.EarlyHandles)
    EARLY_HANDLES = [
        'early_10', 'early_20', 'early_30', 'early_40',
        'early_50', 'early_60', 'early_70', 'early_80', 'early_90',
        'early_95', 'early_97', 'early_100',
    ]

    # Internal names for the 12 Advanced Risk percentile entry widgets
    # (MATLAB: handles.AdvHandles)
    ADV_HANDLES = [
        'adv_10', 'adv_20', 'adv_30', 'adv_40',
        'adv_50', 'adv_60', 'adv_70', 'adv_80', 'adv_90',
        'adv_95', 'adv_97', 'adv_100',
    ]

    # Labels for the 12 percentile positions
    PERCENTILE_LABELS = [
        '10%', '20%', '30%', '40%', '50%', '60%',
        '70%', '80%', '90%', '95%', '97%', '100%',
    ]

    def __init__(self, variables, parent=None):
        ########################################
        ###     Handles Variables            ###
        ########################################

        # Store working copy and backup for cancel/revert
        # MATLAB: handles.Variables    = get(0, 'userdata');
        #         handles.OldVariables = handles.Variables;
        self.variables = variables
        self.old_variables = copy.deepcopy(variables)

        # --- Build the window ------------------------------------------------
        if parent is None:
            self.root = tk.Tk()
        else:
            self.root = tk.Toplevel(parent)

        # MATLAB: set(handles.figure1, 'color', [0.6 0.6 1])
        # MATLAB: set(handles.figure1, 'name', 'Risk graphs', 'NumberTitle','off')
        self.root.title('Risk graphs')
        self.root.configure(bg='#9999ff')  # [0.6 0.6 1] in MATLAB

        # Dictionary to hold all widget references (mirrors MATLAB handles.xxx)
        self.widgets = {}

        self._build_gui()

        # MATLAB: handles = MakeImagesCurrent(hObject, handles);
        self._make_images_current()

    # -----------------------------------------------------------------
    #  GUI construction
    # -----------------------------------------------------------------
    def _build_gui(self):
        """Build all widgets."""
        root = self.root

        # MATLAB Values for percentile positions (used for marker placement)
        # Values = [10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 97, 100]*5
        # These are 1-based MATLAB positions in the 500-element array.

        # ===============================================================
        # Individual Risk section
        # ===============================================================
        ind_frame = tk.LabelFrame(root, text='Individual Risk',
                                  bg='#9999ff')
        ind_frame.pack(padx=5, pady=5, fill='both', expand=True)

        # --- Matplotlib figure for Individual Risk bar chart ---
        self.ind_fig = Figure(figsize=(7, 2.5), dpi=100)
        self.ind_ax = self.ind_fig.add_subplot(111)
        self.ind_canvas = FigureCanvasTkAgg(self.ind_fig, master=ind_frame)
        self.ind_canvas.get_tk_widget().pack(side='top', fill='both', expand=True)

        # --- Percentile entry fields for Individual Risk ---
        ind_entry_frame = tk.Frame(ind_frame, bg='#9999ff')
        ind_entry_frame.pack(side='top', fill='x', padx=2, pady=2)
        for idx, name in enumerate(self.INDIVIDUAL_RISK_HANDLES):
            sub = tk.Frame(ind_entry_frame, bg='#9999ff')
            sub.grid(row=idx // 6, column=idx % 6, padx=2, pady=1)
            tk.Label(sub, text=f'{self.PERCENTILE_LABELS[idx]}:',
                     bg='#9999ff', font=('TkDefaultFont', 8)).pack(side='left')
            entry = tk.Entry(sub, width=8, font=('TkDefaultFont', 8))
            entry.pack(side='left')
            entry.bind('<Return>', lambda e, c=idx: self._risk_percentile_callback(c))
            entry.bind('<FocusOut>', lambda e, c=idx: self._risk_percentile_callback(c))
            self.widgets[name] = entry

        # --- Checkboxes for Individual Risk ---
        ind_check_frame = tk.Frame(ind_frame, bg='#9999ff')
        ind_check_frame.pack(side='top', fill='x', padx=2, pady=2)

        self.log_risk_var = tk.IntVar()
        log_risk_cb = tk.Checkbutton(
            ind_check_frame, text='Log scale',
            variable=self.log_risk_var,
            command=self._log_risk_callback,
            bg='#9999ff', activebackground='#9999ff')
        log_risk_cb.pack(side='left', padx=5)
        self.widgets['log_risk'] = log_risk_cb

        self.normalize_risk_var = tk.IntVar()
        normalize_risk_cb = tk.Checkbutton(
            ind_check_frame, text='Normalize',
            variable=self.normalize_risk_var,
            command=self._normalize_risk_callback,
            bg='#9999ff', activebackground='#9999ff')
        normalize_risk_cb.pack(side='left', padx=5)
        self.widgets['normalize_risk'] = normalize_risk_cb

        # ===============================================================
        # Early Progression Risk section
        # ===============================================================
        early_frame = tk.LabelFrame(root, text='Early Progression Risk',
                                    bg='#9999ff')
        early_frame.pack(padx=5, pady=5, fill='both', expand=True)

        # --- Matplotlib figure for Early Risk bar chart ---
        self.early_fig = Figure(figsize=(7, 2.5), dpi=100)
        self.early_ax = self.early_fig.add_subplot(111)
        self.early_canvas = FigureCanvasTkAgg(self.early_fig, master=early_frame)
        self.early_canvas.get_tk_widget().pack(side='top', fill='both', expand=True)

        # --- Percentile entry fields for Early Risk ---
        early_entry_frame = tk.Frame(early_frame, bg='#9999ff')
        early_entry_frame.pack(side='top', fill='x', padx=2, pady=2)
        for idx, name in enumerate(self.EARLY_HANDLES):
            sub = tk.Frame(early_entry_frame, bg='#9999ff')
            sub.grid(row=idx // 6, column=idx % 6, padx=2, pady=1)
            tk.Label(sub, text=f'{self.PERCENTILE_LABELS[idx]}:',
                     bg='#9999ff', font=('TkDefaultFont', 8)).pack(side='left')
            entry = tk.Entry(sub, width=8, font=('TkDefaultFont', 8))
            entry.pack(side='left')
            entry.bind('<Return>', lambda e, c=idx: self._early_percentile_callback(c))
            entry.bind('<FocusOut>', lambda e, c=idx: self._early_percentile_callback(c))
            self.widgets[name] = entry

        # --- Checkboxes for Early Risk ---
        early_check_frame = tk.Frame(early_frame, bg='#9999ff')
        early_check_frame.pack(side='top', fill='x', padx=2, pady=2)

        self.log_early_var = tk.IntVar()
        log_early_cb = tk.Checkbutton(
            early_check_frame, text='Log scale',
            variable=self.log_early_var,
            command=self._log_early_callback,
            bg='#9999ff', activebackground='#9999ff')
        log_early_cb.pack(side='left', padx=5)
        self.widgets['log_early'] = log_early_cb

        self.normalize_early_var = tk.IntVar()
        normalize_early_cb = tk.Checkbutton(
            early_check_frame, text='Normalize',
            variable=self.normalize_early_var,
            command=self._normalize_early_callback,
            bg='#9999ff', activebackground='#9999ff')
        normalize_early_cb.pack(side='left', padx=5)
        self.widgets['normalize_early'] = normalize_early_cb

        # ===============================================================
        # Advanced Progression Risk section
        # ===============================================================
        adv_frame = tk.LabelFrame(root, text='Advanced Progression Risk',
                                  bg='#9999ff')
        adv_frame.pack(padx=5, pady=5, fill='both', expand=True)

        # --- Matplotlib figure for Advanced Risk bar chart ---
        self.adv_fig = Figure(figsize=(7, 2.5), dpi=100)
        self.adv_ax = self.adv_fig.add_subplot(111)
        self.adv_canvas = FigureCanvasTkAgg(self.adv_fig, master=adv_frame)
        self.adv_canvas.get_tk_widget().pack(side='top', fill='both', expand=True)

        # --- Percentile entry fields for Advanced Risk ---
        adv_entry_frame = tk.Frame(adv_frame, bg='#9999ff')
        adv_entry_frame.pack(side='top', fill='x', padx=2, pady=2)
        for idx, name in enumerate(self.ADV_HANDLES):
            sub = tk.Frame(adv_entry_frame, bg='#9999ff')
            sub.grid(row=idx // 6, column=idx % 6, padx=2, pady=1)
            tk.Label(sub, text=f'{self.PERCENTILE_LABELS[idx]}:',
                     bg='#9999ff', font=('TkDefaultFont', 8)).pack(side='left')
            entry = tk.Entry(sub, width=8, font=('TkDefaultFont', 8))
            entry.pack(side='left')
            entry.bind('<Return>', lambda e, c=idx: self._adv_percentile_callback(c))
            entry.bind('<FocusOut>', lambda e, c=idx: self._adv_percentile_callback(c))
            self.widgets[name] = entry

        # --- Checkboxes for Advanced Risk ---
        adv_check_frame = tk.Frame(adv_frame, bg='#9999ff')
        adv_check_frame.pack(side='top', fill='x', padx=2, pady=2)

        self.log_adv_var = tk.IntVar()
        log_adv_cb = tk.Checkbutton(
            adv_check_frame, text='Log scale',
            variable=self.log_adv_var,
            command=self._log_adv_callback,
            bg='#9999ff', activebackground='#9999ff')
        log_adv_cb.pack(side='left', padx=5)
        self.widgets['log_adv'] = log_adv_cb

        self.normalize_adv_var = tk.IntVar()
        normalize_adv_cb = tk.Checkbutton(
            adv_check_frame, text='Normalize',
            variable=self.normalize_adv_var,
            command=self._normalize_adv_callback,
            bg='#9999ff', activebackground='#9999ff')
        normalize_adv_cb.pack(side='left', padx=5)
        self.widgets['normalize_adv'] = normalize_adv_cb

        # ===============================================================
        # Correlate checkbox and Done button
        # ===============================================================
        bottom_frame = tk.Frame(root, bg='#9999ff')
        bottom_frame.pack(padx=5, pady=5, fill='x')

        self.correlate_var = tk.IntVar()
        correlate_cb = tk.Checkbutton(
            bottom_frame, text='Correlate risks',
            variable=self.correlate_var,
            command=self._correlate_callback,
            bg='#9999ff', activebackground='#9999ff')
        correlate_cb.pack(side='left', padx=5)
        self.widgets['correlate'] = correlate_cb

        tk.Button(bottom_frame, text='Done',
                  command=self._done_callback).pack(side='right', padx=10)

    # -----------------------------------------------------------------
    #  Make Images Current  (update all GUI fields and plots)
    # -----------------------------------------------------------------
    def _make_images_current(self):
        """
        Refresh all GUI widgets and plots from self.variables and recompute
        the interpolated risk arrays via adjust_risk_graph.

        This is the Python equivalent of the MATLAB function
        MakeImagesCurrent(hObject, handles).
        """
        FontSz = 9

        # Recompute the interpolated risk arrays
        # MATLAB: handles.Variables = AdjustRiskGraph(handles.Variables);
        self.variables = adjust_risk_graph(self.variables)

        # MATLAB: Values = [10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 97, 100]*5;
        Values = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 97, 100]) * 5

        # --- normalize checkboxes ---
        # MATLAB: if isequal(handles.Variables.RiskNormalize, 'off'),
        #             set(handles.normalize_risk, 'Value', 0)
        #         else set(handles.normalize_risk, 'Value', 1), end
        if self.variables.get('RiskNormalize', 'off') == 'off':
            self.normalize_risk_var.set(0)
        else:
            self.normalize_risk_var.set(1)

        if self.variables.get('RiskEarlyNormalize', 'off') == 'off':
            self.normalize_early_var.set(0)
        else:
            self.normalize_early_var.set(1)

        if self.variables.get('RiskAdvNormalize', 'off') == 'off':
            self.normalize_adv_var.set(0)
        else:
            self.normalize_adv_var.set(1)

        # ===============================================================
        # Individual Risk plot
        # ===============================================================
        # MATLAB: axes(handles.Individual_risk_axis), cla(gca)
        # MATLAB: bar(handles.Variables.IndividualRisk, 'facecolor', 'b'),
        #         set(gca,'xlim',[1 length(handles.Variables.IndividualRisk)]), hold on
        # MATLAB: line([0 length(...)], [1 1], 'color', 'r')
        ax = self.ind_ax
        ax.clear()
        ind_risk = np.asarray(self.variables['IndividualRisk'])
        # MATLAB bar uses 1-based x positions; in matplotlib we use 0-based
        x_positions = np.arange(len(ind_risk))
        ax.bar(x_positions, ind_risk, color='b', width=1.0)
        ax.set_xlim(0, len(ind_risk))
        ax.axhline(y=1, color='r', linewidth=1)

        # MATLAB: if isequal(handles.Variables.RiskScale, 'linear')
        #             set(gca, 'yscale', 'linear')
        #         else
        #             set(gca, 'yscale', 'log')
        #         end
        if self.variables.get('RiskScale', 'linear') == 'linear':
            ax.set_yscale('linear')
        else:
            ax.set_yscale('log')

        # MATLAB: for x1=1:length(Values)
        #             plot(Values(x1), handles.Variables.Ind_Risk_Percentiles(x1),
        #                  '--rs','LineWidth',1,'MarkerEdgeColor','k',
        #                  'MarkerFaceColor','g','MarkerSize',3)
        #         end
        ind_perc = np.asarray(self.variables['Ind_Risk_Percentiles'], dtype=float)
        for x1 in range(len(Values)):
            # MATLAB Values(x1) is 1-based; matplotlib x is 0-based, so
            # we plot at Values[x1]-1 to match the bar positions
            ax.plot(Values[x1] - 1, ind_perc[x1], 's',
                    markeredgecolor='k', markerfacecolor='g', markersize=3)

        # MATLAB: set(gca, 'xtick', [0 50 100 150 200 250 300 350 400 450 500],
        #     'xticklabel', {'0' '10' '20' '30' '40' '50' '60' '70' '80' '90' '100'})
        ax.set_xticks([0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500])
        ax.set_xticklabels(['0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'])
        ax.tick_params(labelsize=FontSz)
        ax.set_xlabel('percent of patients', fontsize=FontSz)
        ax.set_ylabel('relative risk', fontsize=FontSz)
        self.ind_fig.tight_layout()
        self.ind_canvas.draw()

        # ===============================================================
        # Early Progression Risk plot
        # ===============================================================
        # MATLAB: axes(handles.Early_Polyp_axis), cla(gca)
        ax = self.early_ax
        ax.clear()
        early_risk = np.asarray(self.variables['EarlyRisk'])
        x_positions = np.arange(len(early_risk))
        ax.bar(x_positions, early_risk, color='b', width=1.0)
        ax.set_xlim(0, len(early_risk))
        ax.axhline(y=1, color='r', linewidth=1)

        if self.variables.get('RiskEarlyScale', 'linear') == 'linear':
            ax.set_yscale('linear')
        else:
            ax.set_yscale('log')

        early_perc = np.asarray(self.variables['Early_Risk_Percentiles'], dtype=float)
        for x1 in range(len(Values)):
            ax.plot(Values[x1] - 1, early_perc[x1], 's',
                    markeredgecolor='k', markerfacecolor='g', markersize=3)

        ax.set_xticks([0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500])
        ax.set_xticklabels(['0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'])
        ax.tick_params(labelsize=FontSz)
        ax.set_xlabel('percent of patients', fontsize=FontSz)
        ax.set_ylabel('relative risk', fontsize=FontSz)
        self.early_fig.tight_layout()
        self.early_canvas.draw()

        # ===============================================================
        # Advanced Progression Risk plot
        # ===============================================================
        # MATLAB: axes(handles.Adv_Polyp_axis), cla(gca)
        ax = self.adv_ax
        ax.clear()
        adv_risk = np.asarray(self.variables['AdvRisk'])
        x_positions = np.arange(len(adv_risk))
        ax.bar(x_positions, adv_risk, color='b', width=1.0)
        ax.set_xlim(0, len(adv_risk))
        ax.axhline(y=1, color='r', linewidth=1)

        if self.variables.get('RiskAdvScale', 'linear') == 'linear':
            ax.set_yscale('linear')
        else:
            ax.set_yscale('log')

        adv_perc = np.asarray(self.variables['Adv_Risk_Percentiles'], dtype=float)
        for x1 in range(len(Values)):
            ax.plot(Values[x1] - 1, adv_perc[x1], 's',
                    markeredgecolor='k', markerfacecolor='g', markersize=3)

        ax.set_xticks([0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500])
        ax.set_xticklabels(['0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'])
        ax.tick_params(labelsize=FontSz)
        ax.set_xlabel('percent of patients', fontsize=FontSz)
        ax.set_ylabel('relative risk', fontsize=FontSz)
        self.adv_fig.tight_layout()
        self.adv_canvas.draw()

        # ===============================================================
        # Update all percentile entry fields
        # ===============================================================
        # MATLAB: for x1=1:length(handles.IndividualRiskHandles)
        #             set(handles.(handles.IndividualRiskHandles{x1}), 'string',
        #                 num2str(handles.Variables.Ind_Risk_Percentiles(x1)));
        #         end
        for x1 in range(len(self.INDIVIDUAL_RISK_HANDLES)):
            entry = self.widgets[self.INDIVIDUAL_RISK_HANDLES[x1]]
            entry.delete(0, tk.END)
            entry.insert(0, str(self.variables['Ind_Risk_Percentiles'][x1]))

        for x1 in range(len(self.EARLY_HANDLES)):
            entry = self.widgets[self.EARLY_HANDLES[x1]]
            entry.delete(0, tk.END)
            entry.insert(0, str(self.variables['Early_Risk_Percentiles'][x1]))

        for x1 in range(len(self.ADV_HANDLES)):
            entry = self.widgets[self.ADV_HANDLES[x1]]
            entry.delete(0, tk.END)
            entry.insert(0, str(self.variables['Adv_Risk_Percentiles'][x1]))

        # --- log scale checkboxes ---
        # MATLAB: if isequal(handles.Variables.RiskScale, 'linear'),
        #             set(handles.log_risk, 'Value', 0)
        #         else set(handles.log_risk, 'Value', 1), end
        if self.variables.get('RiskScale', 'linear') == 'linear':
            self.log_risk_var.set(0)
        else:
            self.log_risk_var.set(1)

        if self.variables.get('RiskEarlyScale', 'linear') == 'linear':
            self.log_early_var.set(0)
        else:
            self.log_early_var.set(1)

        if self.variables.get('RiskAdvScale', 'linear') == 'linear':
            self.log_adv_var.set(0)
        else:
            self.log_adv_var.set(1)

        # --- correlate checkbox ---
        # MATLAB: if isequal(handles.Variables.RiskCorrelation, 'on'),
        #             set(handles.correlate, 'Value', 1)
        #         else set(handles.correlate, 'Value', 0), end
        if self.variables.get('RiskCorrelation', 'off') == 'on':
            self.correlate_var.set(1)
        else:
            self.correlate_var.set(0)

    # -----------------------------------------------------------------
    #  Correlate checkbox callback
    # -----------------------------------------------------------------
    def _correlate_callback(self):
        """
        Correlate checkbox callback.

        MATLAB equivalent:
            function correlate_Callback(hObject, eventdata, handles)
            tmp=get(handles.correlate, 'Value');
            if isequal(tmp, 1), handles.Variables.RiskCorrelation = 'on';
            else
                handles.Variables.RiskScale = 'off';
            end
            handles = MakeImagesCurrent(hObject, handles);

        NOTE: The MATLAB code has what appears to be a bug -- when unchecked
        it sets RiskScale = 'off' rather than RiskCorrelation = 'off'.
        We preserve this exact behavior to match the MATLAB code.
        """
        tmp = self.correlate_var.get()
        if tmp == 1:
            self.variables['RiskCorrelation'] = 'on'
        else:
            # MATLAB original sets RiskScale = 'off' (likely a bug, but
            # we preserve the exact behavior)
            self.variables['RiskScale'] = 'off'
        self._make_images_current()

    # -----------------------------------------------------------------
    #  Log scale checkbox callbacks
    # -----------------------------------------------------------------
    def _log_risk_callback(self):
        """
        Log scale checkbox callback for Individual Risk.

        MATLAB equivalent:
            function log_risk_Callback(hObject, eventdata, handles)
            tmp=get(handles.log_risk, 'Value');
            if isequal(tmp, 1), handles.Variables.RiskScale = 'logarithmic';
            else
                handles.Variables.RiskScale = 'linear';
            end
            handles = MakeImagesCurrent(hObject, handles);
        """
        tmp = self.log_risk_var.get()
        if tmp == 1:
            self.variables['RiskScale'] = 'logarithmic'
        else:
            self.variables['RiskScale'] = 'linear'
        self._make_images_current()

    def _log_early_callback(self):
        """
        Log scale checkbox callback for Early Risk.

        MATLAB equivalent:
            function log_early_Callback(hObject, eventdata, handles)
            tmp=get(handles.log_early, 'Value');
            if isequal(tmp, 1), handles.Variables.RiskEarlyScale = 'logarithmic';
            else
                handles.Variables.RiskEarlyScale = 'linear';
            end
            handles = MakeImagesCurrent(hObject, handles);
        """
        tmp = self.log_early_var.get()
        if tmp == 1:
            self.variables['RiskEarlyScale'] = 'logarithmic'
        else:
            self.variables['RiskEarlyScale'] = 'linear'
        self._make_images_current()

    def _log_adv_callback(self):
        """
        Log scale checkbox callback for Advanced Risk.

        MATLAB equivalent:
            function log_adv_Callback(hObject, eventdata, handles)
            tmp=get(handles.log_adv, 'Value');
            if isequal(tmp, 1), handles.Variables.RiskAdvScale = 'logarithmic';
            else
                handles.Variables.RiskAdvScale = 'linear';
            end
            handles = MakeImagesCurrent(hObject, handles);
        """
        tmp = self.log_adv_var.get()
        if tmp == 1:
            self.variables['RiskAdvScale'] = 'logarithmic'
        else:
            self.variables['RiskAdvScale'] = 'linear'
        self._make_images_current()

    # -----------------------------------------------------------------
    #  Normalize checkbox callbacks
    # -----------------------------------------------------------------
    def _normalize_risk_callback(self):
        """
        Normalize checkbox callback for Individual Risk.

        MATLAB equivalent:
            function normalize_risk_Callback(hObject, eventdata, handles)
            tmp=get(handles.normalize_risk, 'Value');
            if isequal(tmp, 1), handles.Variables.RiskNormalize = 'on';
            else
                handles.Variables.RiskNormalize = 'off';
            end
            handles = MakeImagesCurrent(hObject, handles);
        """
        tmp = self.normalize_risk_var.get()
        if tmp == 1:
            self.variables['RiskNormalize'] = 'on'
        else:
            self.variables['RiskNormalize'] = 'off'
        self._make_images_current()

    def _normalize_early_callback(self):
        """
        Normalize checkbox callback for Early Risk.

        MATLAB equivalent:
            function normalize_early_Callback(hObject, eventdata, handles)
            tmp=get(handles.normalize_early, 'Value');
            if isequal(tmp, 1), handles.Variables.RiskEarlyNormalize = 'on';
            else
                handles.Variables.RiskEarlyNormalize = 'off';
            end
            handles = MakeImagesCurrent(hObject, handles);
        """
        tmp = self.normalize_early_var.get()
        if tmp == 1:
            self.variables['RiskEarlyNormalize'] = 'on'
        else:
            self.variables['RiskEarlyNormalize'] = 'off'
        self._make_images_current()

    def _normalize_adv_callback(self):
        """
        Normalize checkbox callback for Advanced Risk.

        MATLAB equivalent:
            function normalize_adv_Callback(hObject, eventdata, handles)
            tmp=get(handles.normalize_adv, 'Value');
            if isequal(tmp, 1), handles.Variables.RiskAdvNormalize = 'on';
            else
                handles.Variables.RiskAdvNormalize = 'off';
            end
            handles = MakeImagesCurrent(hObject, handles);
        """
        tmp = self.normalize_adv_var.get()
        if tmp == 1:
            self.variables['RiskAdvNormalize'] = 'on'
        else:
            self.variables['RiskAdvNormalize'] = 'off'
        self._make_images_current()

    # -----------------------------------------------------------------
    #  Risk percentile callbacks
    # -----------------------------------------------------------------
    def _risk_percentile_callback(self, c):
        """
        Callback for Individual Risk percentile entry at 0-based index c.

        MATLAB equivalent (example for c=0, which is MATLAB c=1):
            c=1; tmp=get(handles.(handles.IndividualRiskHandles{c}), 'string');
            [num, succ] = str2num(tmp);
            if succ, if num >=0, handles.Variables.Ind_Risk_Percentiles(c)=num; end, end
            handles = MakeImagesCurrent(hObject, handles);

        Parameters
        ----------
        c : int
            0-based index into INDIVIDUAL_RISK_HANDLES and
            Ind_Risk_Percentiles (MATLAB c is 1-based).
        """
        tmp = self.widgets[self.INDIVIDUAL_RISK_HANDLES[c]].get()
        num, succ = _safe_str2num(tmp)
        if succ:
            if num >= 0:
                self.variables['Ind_Risk_Percentiles'][c] = num
        self._make_images_current()

    def _early_percentile_callback(self, c):
        """
        Callback for Early Risk percentile entry at 0-based index c.

        MATLAB equivalent (example for c=0, which is MATLAB c=1):
            c=1; tmp=get(handles.(handles.EarlyHandles{c}), 'string');
            [num, succ] = str2num(tmp);
            if succ, if num >=0, handles.Variables.Early_Risk_Percentiles(c)=num; end, end
            handles = MakeImagesCurrent(hObject, handles);

        Parameters
        ----------
        c : int
            0-based index into EARLY_HANDLES and
            Early_Risk_Percentiles (MATLAB c is 1-based).
        """
        tmp = self.widgets[self.EARLY_HANDLES[c]].get()
        num, succ = _safe_str2num(tmp)
        if succ:
            if num >= 0:
                self.variables['Early_Risk_Percentiles'][c] = num
        self._make_images_current()

    def _adv_percentile_callback(self, c):
        """
        Callback for Advanced Risk percentile entry at 0-based index c.

        MATLAB equivalent (example for c=0, which is MATLAB c=1):
            c=1; tmp=get(handles.(handles.AdvHandles{c}), 'string');
            [num, succ] = str2num(tmp);
            if succ, if num >=0, handles.Variables.Adv_Risk_Percentiles(c)=num; end, end
            handles = MakeImagesCurrent(hObject, handles);

        Parameters
        ----------
        c : int
            0-based index into ADV_HANDLES and
            Adv_Risk_Percentiles (MATLAB c is 1-based).
        """
        tmp = self.widgets[self.ADV_HANDLES[c]].get()
        num, succ = _safe_str2num(tmp)
        if succ:
            if num >= 0:
                self.variables['Adv_Risk_Percentiles'][c] = num
        self._make_images_current()

    # -----------------------------------------------------------------
    #  Done callback
    # -----------------------------------------------------------------
    def _done_callback(self):
        """
        Done button callback.

        MATLAB equivalent:
            function done_Callback(hObject, eventdata, handles)
            Answer = questdlg('Do you want to keep the settings?', ...
                              'Return?', 'Yes', 'No', 'Cancel', 'Yes');
            if isequal(Answer, 'Cancel')
                return
            elseif isequal(Answer, 'No')
                handles.Variables = handles.OldVariables;
            end
            set(0, 'userdata', handles.Variables);
            uiresume(handles.figure1);
            if ishandle(handles.figure1)
                delete(handles.figure1);
            end
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

        # MATLAB: set(0, 'userdata', handles.Variables);
        # In Python, self.variables is already the shared reference.

        # MATLAB: uiresume(handles.figure1);
        # MATLAB: if ishandle(handles.figure1), delete(handles.figure1); end
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

def risk_settings(variables):
    """
    Open the Risk Settings dialog and return the (possibly modified)
    Variables dictionary.

    This is the main entry point, equivalent to calling the MATLAB function
    ``Risk_Settings`` which opens a GUIDE dialog, blocks until the user
    clicks Done, and returns the updated userdata.

    Parameters
    ----------
    variables : dict
        The simulation Variables dictionary.  Must contain the following keys:
        - 'Ind_Risk_Percentiles' : list of 12 floats
        - 'Early_Risk_Percentiles' : list of 12 floats
        - 'Adv_Risk_Percentiles' : list of 12 floats
        - 'IndividualRisk' : list of 500 floats (recomputed by adjust_risk_graph)
        - 'EarlyRisk' : list of 500 floats (recomputed by adjust_risk_graph)
        - 'AdvRisk' : list of 500 floats (recomputed by adjust_risk_graph)
        - 'RiskNormalize' : 'on' or 'off'
        - 'RiskEarlyNormalize' : 'on' or 'off'
        - 'RiskAdvNormalize' : 'on' or 'off'
        - 'RiskScale' : 'linear' or 'logarithmic'
        - 'RiskEarlyScale' : 'linear' or 'logarithmic'
        - 'RiskAdvScale' : 'linear' or 'logarithmic'
        - 'RiskCorrelation' : 'on' or 'off'
        Modified in-place if the user accepts changes.

    Returns
    -------
    dict
        The same ``variables`` dict reference (possibly reverted if the user
        chose 'No').
    """
    gui = RiskSettingsGUI(variables)
    gui.run()
    return variables


# ===================================================================
#  STANDALONE EXECUTION
# ===================================================================

if __name__ == '__main__':
    # Example: create a minimal Variables dict for testing.
    # This mirrors the risk-related structure expected by the GUI,
    # using values from CMOST19 settings.
    test_variables = {
        'Ind_Risk_Percentiles': [0.0, 0.0, 0.05, 0.7, 0.9, 1.4,
                                  1.45, 1.52, 1.58, 1.8, 7.0, 12.0],
        'Early_Risk_Percentiles': [0.05, 0.2, 0.2, 0.3, 0.5, 1.0,
                                    1.1, 1.6, 1.8, 2.0, 2.7, 4.0],
        'Adv_Risk_Percentiles': [0.0, 0.0, 0.0, 0.2, 0.25, 0.3,
                                  0.5, 0.7, 0.9, 1.2, 1.7, 2.5],
        'IndividualRisk': [0.0] * 500,
        'EarlyRisk': [0.0] * 500,
        'AdvRisk': [0.0] * 500,
        'RiskNormalize': 'off',
        'RiskEarlyNormalize': 'off',
        'RiskAdvNormalize': 'off',
        'RiskScale': 'linear',
        'RiskEarlyScale': 'linear',
        'RiskAdvScale': 'linear',
        'RiskCorrelation': 'on',
    }

    result = risk_settings(test_variables)
    print('Final Ind_Risk_Percentiles:', result['Ind_Risk_Percentiles'])
    print('Final Early_Risk_Percentiles:', result['Early_Risk_Percentiles'])
    print('Final Adv_Risk_Percentiles:', result['Adv_Risk_Percentiles'])
    print('IndividualRisk (first 10):', result['IndividualRisk'][:10])
    print('EarlyRisk (first 10):', result['EarlyRisk'][:10])
    print('AdvRisk (first 10):', result['AdvRisk'][:10])
