import tkinter as tk
window = tk.Tk()

greeting = tk.Label(text="Hello, Tkinter")

run=tk.Label(text='run a molpro calculation')
run.grid(column=0, row=0)
runentry = tk.Entry(fg="yellow", bg="yellow", width=10)
runentry.grid(column=1, row=0)

window.mainloop()