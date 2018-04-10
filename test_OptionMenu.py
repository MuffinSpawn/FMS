# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 15:03:15 2018

@author: plane
"""

import tkinter as tk
import tkinter.ttk as ttk

def let_it_out(selections):
    for selection in selections:
        print(selection.get())
        print(selection.get())

def main():
    root = tk.Tk()

    notebook = ttk.Notebook(root)
    notebook.grid(row=0,column=1, rowspan=6, columnspan=3)

    page = tk.Frame(notebook)
    page.grid(row=0, column=0)
    notebook.add(page, text="A Tab")

    choices = ['A', 'B', 'C']
    selections = []
    selection = tk.StringVar()
    selection.set(choices[0])
    selections.append(selection)
    tk.OptionMenu(page, selection, *choices).grid(row=0, column=0)
    tk.Button(page, text='Let It Out!', command=lambda: let_it_out(selections)).grid(row=0, column=1)

    root.mainloop()

if __name__ == '__main__':
    main()