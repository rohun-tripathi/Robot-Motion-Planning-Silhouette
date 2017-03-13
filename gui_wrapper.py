#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

import ScrolledText
import Tkinter

import main

SCREEN_BACKGROUND = "dodgerblue3"
LABEL_BACKGROUND = "darkorchid"

SPACE = " "

SET_DEFAULT_DIMENSIONS = 3
SET_DEFAULT_NUMBER = 4
ALLOWED_DIMENSION_VALUES = [3]


class path_finder_tk(Tkinter.Tk):
    def __init__(self, parent):
        Tkinter.Tk.__init__(self, parent)
        self.parent = parent
        self.row_number = -1
        self.initialize()

    def increment_row_and_add_field_label(self, text, pady=8):
        label = Tkinter.Label(self, text=text, anchor="w", fg="white", bg=LABEL_BACKGROUND)
        label.grid(column=0, row=self.increment_row_number(), columnspan=1, sticky='EW', padx=5, pady=pady)

    def add_major_label(self, txt):
        headinglabel = Tkinter.Label(self, text=txt, anchor="center", fg="white", bg=LABEL_BACKGROUND)
        headinglabel.grid(column=0, row=self.increment_row_number(), sticky='EW', columnspan=4, padx=4, pady=16)

    def increment_row_number(self):
        self.row_number += 1
        return self.row_number

    def set_console_label_value(self, value):
        self.label.config(state='normal')
        self.label.delete("0.0", "end")
        self.label.insert("insert", value)
        self.label.config(state='disabled')

    def add_ellipse_to_input(self, dimension, input_lines, value):
        lines = value.strip("\n").split("\n")
        assert len(lines) == 1 + dimension + 1

        concatenated_line = ""
        for line in lines:
            concatenated_line += line + SPACE
        input_lines.append(concatenated_line.strip())

    def initialize(self):
        self.grid()

        self.add_major_label(u"Path Finder")

        self.increment_row_and_add_field_label(u"Number of Dimensions : ")
        self.dimEntryVariable = Tkinter.StringVar()
        self.dimEntry = Tkinter.Entry(self, textvariable=self.dimEntryVariable)
        self.dimEntry.grid(column=1, row=self.row_number, sticky='EW', padx=5, pady=8)
        self.dimEntryVariable.set(SET_DEFAULT_DIMENSIONS)

        self.increment_row_and_add_field_label(u"Bounding Arena : ")
        self.primary_ellipsoid_entry = ScrolledText.ScrolledText(self, wrap="word", height=5)
        self.primary_ellipsoid_entry.grid(column=1, row=self.row_number, sticky='W', padx=5, pady=8)
        self.primary_ellipsoid_entry.insert('insert', '''0 0 0\n1 0 0\n0 1 0\n0 0 1\n20 20 20''')

        self.increment_row_and_add_field_label(u"Obstacles : ")
        self.obstacles_entry = ScrolledText.ScrolledText(self, wrap="word", height=10)
        self.obstacles_entry.grid(column=1, row=self.row_number, sticky='W', padx=5, pady=8)
        self.obstacles_entry.insert('insert',
                                    '''-4 0 0\n1 0 0\n0 1 0\n0 0 1\n3 3 3\n\n4 0 0\n1 0 0\n0 1 0\n0 0 1\n3 3 3''')

        self.increment_row_and_add_field_label(u"Save to File : ")
        self.file_name_variable = Tkinter.StringVar()
        self.file_name_entry = Tkinter.Entry(self, textvariable=self.file_name_variable)
        self.file_name_entry.grid(column=1, row=self.row_number, sticky='EW', padx=5, pady=8)
        self.file_name_variable.set("/home/riot/Desktop/roadmap")

        button = Tkinter.Button(self, text=u"Draw Silhouette Curves", command=self.on_compute_road_map)
        button.grid(column=0, row=self.increment_row_number())

        self.increment_row_and_add_field_label(u"Console : ")
        self.label = ScrolledText.ScrolledText(self, wrap="word", fg="black", bg="white", height=5)
        self.label.grid(column=1, row=self.row_number, columnspan=4, sticky='EW', padx=5, pady=8)

        self.grid_columnconfigure(0, weight=1)
        self.resizable(True, False)
        self.update()
        self.geometry(self.geometry())
        self.configure(background=SCREEN_BACKGROUND)

        self.dimEntry.focus_set()
        self.dimEntry.selection_range(0, Tkinter.END)

    def on_compute_road_map(self):
        input_lines = []
        self.set_console_label_value("")

        try:
            dimension = int(self.dimEntryVariable.get())
            assert dimension in ALLOWED_DIMENSION_VALUES
            input_lines.append(str(dimension))
        except AssertionError:
            self.set_console_label_value("Allowed Dimension values are : " + str(ALLOWED_DIMENSION_VALUES))
            self.dimEntry.focus_set()
            self.dimEntry.selection_range(0, Tkinter.END)
            return
        except Exception:
            self.set_console_label_value("Unknown Error in Dimension value. Please check the parameters")
            self.dimEntry.focus_set()
            self.dimEntry.selection_range(0, Tkinter.END)
            return

        try:
            value = self.primary_ellipsoid_entry.get("0.0", Tkinter.END)
            self.add_ellipse_to_input(dimension, input_lines, value)
        except AssertionError:
            self.set_console_label_value("Primary Ellipsoid Input : Wrong number of lines. Should be 1 + dimension + 1")
            self.primary_ellipsoid_entry.focus_set()
            self.primary_ellipsoid_entry.selection_get(0)
            return
        except Exception:
            self.set_console_label_value("Unknown Error in Primary Ellipsoid value. Please check the parameters")
            self.primary_ellipsoid_entry.focus_set()
            self.primary_ellipsoid_entry.selection_get(0)
            return

        try:
            value = self.obstacles_entry.get("0.0", Tkinter.END)
            obstacle_list = value.split("\n\n")
            for obstacle in obstacle_list:
                self.add_ellipse_to_input(dimension, input_lines, obstacle)
        except AssertionError:
            self.set_console_label_value("Obstacles Input : Wrong number of lines. Should be 1 + dimension + 1")
            self.obstacles_entry.focus_set()
            self.obstacles_entry.selection_get(0)
            return
        except Exception:
            self.set_console_label_value("Unknown Error in Obstacles input. Please check the parameters")
            self.obstacles_entry.focus_set()
            self.obstacles_entry.selection_get(0)
            return

        try:
            file_name = str(self.file_name_entry.get()).strip()
            file_name = "foo" if file_name == "" else file_name

            main.run_silhouette_method(input_lines, file_name)
            self.set_console_label_value("Road Map computation Done.")
        except Exception as e:
            self.set_console_label_value("Could Not Compute Road Map. Please check parameter values.")

        self.dimEntry.focus_set()
        self.dimEntry.selection_range(0, Tkinter.END)


if __name__ == "__main__":
    app = path_finder_tk(None)
    app.title("Silhouette Method")
    app.mainloop()