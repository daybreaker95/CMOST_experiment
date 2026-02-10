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
# MATLAB's  handles.Variables  is stored as  self.variables  (a dict).
# MATLAB's  handles.Variables.Cost  becomes  self.variables['Cost']  (a dict).
# MATLAB structs become Python dicts, consistent with NumberCrunching_100000.py,
# Evaluation.py, and Colonoscopy_settings.py.
#
# The Cost_Settings GUI edits 29 cost fields plus a "near future" checkbox.
# When NearFuture is 'off', the current cost variable names are used.
# When NearFuture is 'on', the future cost variable names (prefixed Fut)
# are used for the cancer treatment cost entries (indices 8..23).
#
# Fieldhandles (GUI widget names) map 1:1 to VariableHandles (Cost dict keys).
# The mapping depends on the NearFuture toggle state.
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
#  COST SETTINGS GUI
# ===================================================================

class CostSettingsGUI:
    """
    Python/tkinter equivalent of the MATLAB GUIDE Cost_Settings dialog.

    The GUI allows the user to view and edit cost parameters stored in
    Variables['Cost'].  There are 29 editable fields covering:
      - Colonoscopy and colonoscopy with polyp removal costs
      - Rectosigmoidoscopy and rectosigmoidoscopy with polyp removal costs
      - Complication costs (perforation, serosa burn, bleeding, severe bleeding)
      - Cancer treatment costs by stage (I-IV) for initial, continuing,
        terminal, and terminal opportunity costs
      - Screening test costs (FOBT, I-FOBT, Sept9 variants, other)

    A "Near Future" checkbox toggles between current and future cost
    variable names for the cancer treatment entries.

    Changes are applied to the Variables dict.  The user may accept,
    reject, or cancel changes via the Done button.

    Parameters
    ----------
    variables : dict
        The simulation Variables dictionary.  Must contain a 'Cost'
        sub-dictionary.  This dictionary is modified in-place if the
        user accepts changes.
    parent : tk.Tk or tk.Toplevel, optional
        Parent window.  If None, a new Tk root is created.
    """

    # these variables are internal references for the elements in the GUI
    # (MATLAB: handles.Fieldhandles)
    # These are the GUI widget names - 29 fields total
    FIELD_HANDLES = [
        'Colonoscopy', 'Colonoscopy_Polyp', 'Rectosigmoidoscopy', 'Rectosigmoidoscopy_Polyp',
        'Perforation', 'serosa_burn', 'bleeding', 'severe_bleeding_transfusion',
        'Cancer_I_initial', 'Cancer_II_initial', 'Cancer_III_initial', 'Cancer_IV_initial',
        'cancer_I_cont', 'Cancer_II_cont', 'Cancer_III_cont', 'Cancer_IV_cont',
        'Cancer_I_terminal', 'Cancer_II_terminal', 'Cancer_III_terminal', 'Cancer_VI_terminal',
        'Cancer_I_terminal_oc', 'Cancer_II_terminal_oc', 'Cancer_III_terminal_oc', 'Cancer_IV_terminal_oc',
        'FOBT', 'I_FOBT', 'Sept9_HighSens', 'Sept9_HighSpec', 'other',
    ]

    # Variable handles when NearFuture is 'off' (current costs)
    # (MATLAB: handles.VariableHandles in the 'off' branch)
    VARIABLE_HANDLES_CURRENT = [
        'Colonoscopy', 'Colonoscopy_Polyp', 'Sigmoidoscopy', 'Sigmoidoscopy_Polyp',
        'Colonoscopy_Perforation', 'Colonoscopy_Serosal_burn', 'Colonoscopy_bleed', 'Colonoscopy_bleed_transfusion',
        'Initial_I', 'Initial_II', 'Initial_III', 'Initial_IV',
        'Cont_I', 'Cont_II', 'Cont_III', 'Cont_IV',
        'Final_I', 'Final_II', 'Final_III', 'Final_IV',
        'Final_oc_I', 'Final_oc_II', 'Final_oc_III', 'Final_oc_IV',
        'FOBT', 'I_FOBT', 'Sept9_HighSens', 'Sept9_HighSpec', 'other',
    ]

    # Variable handles when NearFuture is 'on' (future costs)
    # (MATLAB: handles.VariableHandles in the 'on' branch)
    VARIABLE_HANDLES_FUTURE = [
        'Colonoscopy', 'Colonoscopy_Polyp', 'Sigmoidoscopy', 'Sigmoidoscopy_Polyp',
        'Colonoscopy_Perforation', 'Colonoscopy_Serosal_burn', 'Colonoscopy_bleed', 'Colonoscopy_bleed_transfusion',
        'FutInitial_I', 'FutInitial_II', 'FutInitial_III', 'FutInitial_IV',
        'FutCont_I', 'FutCont_II', 'FutCont_III', 'FutCont_IV',
        'FutFinal_I', 'FutFinal_II', 'FutFinal_III', 'FutFinal_IV',
        'FutFinal_oc_I', 'FutFinal_oc_II', 'FutFinal_oc_III', 'FutFinal_oc_IV',
        'FOBT', 'I_FOBT', 'Sept9_HighSens', 'Sept9_HighSpec', 'other',
    ]

    # Display labels for the GUI (human-readable names for each field)
    FIELD_LABELS = [
        'Colonoscopy', 'Colonoscopy + Polyp', 'Rectosigmoidoscopy', 'Rectosigmoidoscopy + Polyp',
        'Perforation', 'Serosa burn', 'Bleeding', 'Severe bleeding / transfusion',
        'Cancer I initial', 'Cancer II initial', 'Cancer III initial', 'Cancer IV initial',
        'Cancer I continuing', 'Cancer II continuing', 'Cancer III continuing', 'Cancer IV continuing',
        'Cancer I terminal', 'Cancer II terminal', 'Cancer III terminal', 'Cancer IV terminal',
        'Cancer I terminal (opp. cost)', 'Cancer II terminal (opp. cost)', 'Cancer III terminal (opp. cost)', 'Cancer IV terminal (opp. cost)',
        'FOBT', 'I-FOBT', 'Sept9 High Sens', 'Sept9 High Spec', 'Other',
    ]

    def __init__(self, variables, parent=None):
        ########################################
        ###     Handles Variables            ###
        ########################################

        # Store working copy and backup for cancel/revert
        # MATLAB: handles.Variables = get(0, 'userdata');
        #         handles.OldVariables = handles.Variables;
        self.variables = variables
        self.old_variables = copy.deepcopy(variables)

        # --- Build the window ------------------------------------------------
        if parent is None:
            self.root = tk.Tk()
        else:
            self.root = tk.Toplevel(parent)

        # MATLAB: set(handles.figure1, 'color', [0.6 0.6 1])
        # MATLAB: set(handles.figure1, 'name', 'Costs', 'NumberTitle','off')
        self.root.title('Costs')
        self.root.configure(bg='#9999ff')  # [0.6 0.6 1] in MATLAB

        # Dictionary to hold all widget references (mirrors MATLAB handles.xxx)
        self.widgets = {}

        self._build_gui()

        # MATLAB: handles = MakeImagesCurrent(hObject, handles);
        self._make_images_current()

    # -----------------------------------------------------------------
    #  Get current variable handles based on NearFuture state
    # -----------------------------------------------------------------
    def _get_variable_handles(self):
        """
        Return the appropriate variable handle list based on NearFuture state.

        MATLAB equivalent:
            if isequal(handles.Variables.Cost.NearFuture, 'off')
                handles.VariableHandles = { ... current names ... };
            else
                handles.VariableHandles = { ... future names ... };
            end
        """
        if self.variables['Cost'].get('NearFuture', 'off') == 'off':
            return self.VARIABLE_HANDLES_CURRENT
        else:
            return self.VARIABLE_HANDLES_FUTURE

    # -----------------------------------------------------------------
    #  GUI construction
    # -----------------------------------------------------------------
    def _build_gui(self):
        """Build all widgets."""
        root = self.root

        # --- Near Future checkbox --------------------------------------------
        checkbox_frame = tk.Frame(root, bg='#9999ff')
        checkbox_frame.pack(padx=5, pady=5, fill='x')

        self.near_future_var = tk.IntVar()
        # Initialize checkbox from current NearFuture state
        if self.variables['Cost'].get('NearFuture', 'off') == 'on':
            self.near_future_var.set(1)
        else:
            self.near_future_var.set(0)

        near_future_cb = tk.Checkbutton(
            checkbox_frame, text='Near Future Costs',
            variable=self.near_future_var,
            command=self._near_future_callback,
            bg='#9999ff', activebackground='#9999ff')
        near_future_cb.pack(side='left')
        self.widgets['near_future'] = near_future_cb

        # --- Procedure Costs section (indices 0-3) ---------------------------
        proc_frame = tk.LabelFrame(root, text='Procedure Costs', bg='#9999ff')
        proc_frame.pack(padx=5, pady=5, fill='x')

        for idx in range(4):
            row = idx // 2
            col = idx % 2
            sub = tk.Frame(proc_frame, bg='#9999ff')
            sub.grid(row=row, column=col, padx=4, pady=2, sticky='w')
            tk.Label(sub, text=f'{self.FIELD_LABELS[idx]}:', bg='#9999ff').pack(side='left')
            entry = tk.Entry(sub, width=12)
            entry.pack(side='left', padx=2)
            # Bind callback on Return key and FocusOut
            entry.bind('<Return>', lambda e, c=idx: self._field_callback(c))
            entry.bind('<FocusOut>', lambda e, c=idx: self._field_callback(c))
            self.widgets[self.FIELD_HANDLES[idx]] = entry

        # --- Complication Costs section (indices 4-7) ------------------------
        comp_frame = tk.LabelFrame(root, text='Complication Costs', bg='#9999ff')
        comp_frame.pack(padx=5, pady=5, fill='x')

        for idx in range(4, 8):
            row = (idx - 4) // 2
            col = (idx - 4) % 2
            sub = tk.Frame(comp_frame, bg='#9999ff')
            sub.grid(row=row, column=col, padx=4, pady=2, sticky='w')
            tk.Label(sub, text=f'{self.FIELD_LABELS[idx]}:', bg='#9999ff').pack(side='left')
            entry = tk.Entry(sub, width=12)
            entry.pack(side='left', padx=2)
            entry.bind('<Return>', lambda e, c=idx: self._field_callback(c))
            entry.bind('<FocusOut>', lambda e, c=idx: self._field_callback(c))
            self.widgets[self.FIELD_HANDLES[idx]] = entry

        # --- Cancer Initial Treatment Costs section (indices 8-11) -----------
        initial_frame = tk.LabelFrame(root, text='Cancer Initial Treatment Costs',
                                       bg='#9999ff')
        initial_frame.pack(padx=5, pady=5, fill='x')

        for idx in range(8, 12):
            col = idx - 8
            sub = tk.Frame(initial_frame, bg='#9999ff')
            sub.grid(row=0, column=col, padx=4, pady=2, sticky='w')
            tk.Label(sub, text=f'{self.FIELD_LABELS[idx]}:', bg='#9999ff').pack(side='left')
            entry = tk.Entry(sub, width=10)
            entry.pack(side='left', padx=2)
            entry.bind('<Return>', lambda e, c=idx: self._field_callback(c))
            entry.bind('<FocusOut>', lambda e, c=idx: self._field_callback(c))
            self.widgets[self.FIELD_HANDLES[idx]] = entry

        # --- Cancer Continuing Treatment Costs section (indices 12-15) -------
        cont_frame = tk.LabelFrame(root, text='Cancer Continuing Treatment Costs',
                                    bg='#9999ff')
        cont_frame.pack(padx=5, pady=5, fill='x')

        for idx in range(12, 16):
            col = idx - 12
            sub = tk.Frame(cont_frame, bg='#9999ff')
            sub.grid(row=0, column=col, padx=4, pady=2, sticky='w')
            tk.Label(sub, text=f'{self.FIELD_LABELS[idx]}:', bg='#9999ff').pack(side='left')
            entry = tk.Entry(sub, width=10)
            entry.pack(side='left', padx=2)
            entry.bind('<Return>', lambda e, c=idx: self._field_callback(c))
            entry.bind('<FocusOut>', lambda e, c=idx: self._field_callback(c))
            self.widgets[self.FIELD_HANDLES[idx]] = entry

        # --- Cancer Terminal Treatment Costs section (indices 16-19) ---------
        term_frame = tk.LabelFrame(root, text='Cancer Terminal Treatment Costs',
                                    bg='#9999ff')
        term_frame.pack(padx=5, pady=5, fill='x')

        for idx in range(16, 20):
            col = idx - 16
            sub = tk.Frame(term_frame, bg='#9999ff')
            sub.grid(row=0, column=col, padx=4, pady=2, sticky='w')
            tk.Label(sub, text=f'{self.FIELD_LABELS[idx]}:', bg='#9999ff').pack(side='left')
            entry = tk.Entry(sub, width=10)
            entry.pack(side='left', padx=2)
            entry.bind('<Return>', lambda e, c=idx: self._field_callback(c))
            entry.bind('<FocusOut>', lambda e, c=idx: self._field_callback(c))
            self.widgets[self.FIELD_HANDLES[idx]] = entry

        # --- Cancer Terminal Opportunity Costs section (indices 20-23) -------
        term_oc_frame = tk.LabelFrame(root, text='Cancer Terminal Opportunity Costs',
                                       bg='#9999ff')
        term_oc_frame.pack(padx=5, pady=5, fill='x')

        for idx in range(20, 24):
            col = idx - 20
            sub = tk.Frame(term_oc_frame, bg='#9999ff')
            sub.grid(row=0, column=col, padx=4, pady=2, sticky='w')
            tk.Label(sub, text=f'{self.FIELD_LABELS[idx]}:', bg='#9999ff').pack(side='left')
            entry = tk.Entry(sub, width=10)
            entry.pack(side='left', padx=2)
            entry.bind('<Return>', lambda e, c=idx: self._field_callback(c))
            entry.bind('<FocusOut>', lambda e, c=idx: self._field_callback(c))
            self.widgets[self.FIELD_HANDLES[idx]] = entry

        # --- Screening Test Costs section (indices 24-28) --------------------
        screen_frame = tk.LabelFrame(root, text='Screening Test Costs',
                                      bg='#9999ff')
        screen_frame.pack(padx=5, pady=5, fill='x')

        for idx in range(24, 29):
            col = idx - 24
            sub = tk.Frame(screen_frame, bg='#9999ff')
            sub.grid(row=0, column=col, padx=4, pady=2, sticky='w')
            tk.Label(sub, text=f'{self.FIELD_LABELS[idx]}:', bg='#9999ff').pack(side='left')
            entry = tk.Entry(sub, width=10)
            entry.pack(side='left', padx=2)
            entry.bind('<Return>', lambda e, c=idx: self._field_callback(c))
            entry.bind('<FocusOut>', lambda e, c=idx: self._field_callback(c))
            self.widgets[self.FIELD_HANDLES[idx]] = entry

        # --- Done button -----------------------------------------------------
        btn_frame = tk.Frame(root, bg='#9999ff')
        btn_frame.pack(padx=5, pady=10)
        tk.Button(btn_frame, text='Done', command=self._done_callback).pack()

    # -----------------------------------------------------------------
    #  Make Images Current  (update all GUI fields from Variables)
    # -----------------------------------------------------------------
    def _make_images_current(self):
        """
        Refresh all GUI widgets from self.variables['Cost'] based on
        the current NearFuture state.

        This is the Python equivalent of the MATLAB function
        MakeImagesCurrent(hObject, handles).

        MATLAB equivalent:
            if isequal(handles.Variables.Cost.NearFuture, 'off')
                handles.VariableHandles = { ... current ... };
            else
                handles.VariableHandles = { ... future ... };
            end
            for f=1:length(handles.Fieldhandles)
                set(handles.(handles.Fieldhandles{f}), 'string',
                    num2str(handles.Variables.Cost.(handles.VariableHandles{f})));
            end
        """
        variable_handles = self._get_variable_handles()
        cost = self.variables['Cost']

        # adjust graphs for new polyps
        # MATLAB: for f=1:length(handles.Fieldhandles)
        #             set(handles.(handles.Fieldhandles{f}), 'string',
        #                 num2str(handles.Variables.Cost.(handles.VariableHandles{f})));
        #         end
        for f in range(len(self.FIELD_HANDLES)):
            field_name = self.FIELD_HANDLES[f]
            var_name = variable_handles[f]
            entry = self.widgets[field_name]
            entry.delete(0, tk.END)
            entry.insert(0, str(cost.get(var_name, 0)))

    # -----------------------------------------------------------------
    #  Field callbacks  (generic for all 29 cost fields)
    # -----------------------------------------------------------------
    def _field_callback(self, c):
        """
        Generic callback for cost field at 0-based index c.

        This is the Python equivalent of each individual MATLAB callback:
            function Colonoscopy_Callback(hObject, eventdata, handles)
                c=1; tmp=get(handles.(handles.Fieldhandles{c}), 'string');
                [num, succ] = str2num(tmp);
                if succ, if num >=0,
                    handles.Variables.Cost.(handles.VariableHandles{c})=num;
                end, end,
                MakeImagesCurrent(hObject, handles);

        All 29 callbacks in the MATLAB code follow exactly the same pattern,
        differing only in the value of c (1..29 in MATLAB, 0..28 in Python).

        Parameters
        ----------
        c : int
            0-based index into FIELD_HANDLES / VariableHandles.
        """
        variable_handles = self._get_variable_handles()
        field_name = self.FIELD_HANDLES[c]
        var_name = variable_handles[c]

        tmp = self.widgets[field_name].get()
        num, succ = _safe_str2num(tmp)
        if succ:
            if num >= 0:
                self.variables['Cost'][var_name] = num
        self._make_images_current()

    # -----------------------------------------------------------------
    #  Near Future checkbox callback
    # -----------------------------------------------------------------
    def _near_future_callback(self):
        """
        Near Future checkbox callback.

        MATLAB equivalent:
            function near_future_Callback(hObject, eventdata, handles)
            tmp=get(handles.near_future, 'Value');
            if isequal(tmp, 1), handles.Variables.Cost.NearFuture = 'on';
            else
                 handles.Variables.Cost.NearFuture = 'off';
            end
            handles = MakeImagesCurrent(hObject, handles);
        """
        tmp = self.near_future_var.get()
        if tmp == 1:
            self.variables['Cost']['NearFuture'] = 'on'
        else:
            self.variables['Cost']['NearFuture'] = 'off'
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

def cost_settings(variables):
    """
    Open the Cost Settings dialog and return the (possibly modified)
    Variables dictionary.

    This is the main entry point, equivalent to calling the MATLAB function
    ``Cost_Settings`` which opens a GUIDE dialog, blocks until the
    user clicks Done, and returns the updated userdata.

    Parameters
    ----------
    variables : dict
        The simulation Variables dictionary.  Must contain a 'Cost'
        sub-dictionary.  Modified in-place if the user accepts changes.

    Returns
    -------
    dict
        The same ``variables`` dict reference (possibly reverted if the user
        chose 'No').
    """
    gui = CostSettingsGUI(variables)
    gui.run()
    return variables


# ===================================================================
#  STANDALONE EXECUTION
# ===================================================================

if __name__ == '__main__':
    # Example: create a minimal Variables dict for testing
    # This mirrors the Cost sub-dictionary structure expected by the GUI.
    test_cost = {
        'NearFuture': 'off',
        # Procedure costs
        'Colonoscopy': 1563,
        'Colonoscopy_Polyp': 2010,
        'Sigmoidoscopy': 478,
        'Sigmoidoscopy_Polyp': 925,
        # Complication costs
        'Colonoscopy_Perforation': 11474,
        'Colonoscopy_Serosal_burn': 5000,
        'Colonoscopy_bleed': 3312,
        'Colonoscopy_bleed_transfusion': 8000,
        # Cancer initial treatment costs (current)
        'Initial_I': 32392,
        'Initial_II': 40002,
        'Initial_III': 52424,
        'Initial_IV': 72749,
        # Cancer continuing treatment costs (current)
        'Cont_I': 2177,
        'Cont_II': 3121,
        'Cont_III': 5498,
        'Cont_IV': 25710,
        # Cancer terminal treatment costs (current)
        'Final_I': 36046,
        'Final_II': 44984,
        'Final_III': 57482,
        'Final_IV': 60281,
        # Cancer terminal opportunity costs (current)
        'Final_oc_I': 36046,
        'Final_oc_II': 44984,
        'Final_oc_III': 57482,
        'Final_oc_IV': 60281,
        # Cancer initial treatment costs (future)
        'FutInitial_I': 32392,
        'FutInitial_II': 40002,
        'FutInitial_III': 52424,
        'FutInitial_IV': 72749,
        # Cancer continuing treatment costs (future)
        'FutCont_I': 2177,
        'FutCont_II': 3121,
        'FutCont_III': 5498,
        'FutCont_IV': 25710,
        # Cancer terminal treatment costs (future)
        'FutFinal_I': 36046,
        'FutFinal_II': 44984,
        'FutFinal_III': 57482,
        'FutFinal_IV': 60281,
        # Cancer terminal opportunity costs (future)
        'FutFinal_oc_I': 36046,
        'FutFinal_oc_II': 44984,
        'FutFinal_oc_III': 57482,
        'FutFinal_oc_IV': 60281,
        # Screening test costs
        'FOBT': 22,
        'I_FOBT': 30,
        'Sept9_HighSens': 150,
        'Sept9_HighSpec': 150,
        'other': 0,
    }
    test_variables = {'Cost': test_cost}

    result = cost_settings(test_variables)
    print('Final Cost settings:')
    for key, val in result['Cost'].items():
        print(f'  {key}: {val}')
