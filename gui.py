"""Local graphical entry point; computation is shared with the CLI."""
import queue
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from app_version import VERSION
from network_core import analyze, EXPLORATION_NOTICE, POOLED_NOTICE


def build_options(counts, output, metadata, sample_column, group_column, prevalence, abundance, rho, permutations, repeats, independent, taxonomy='', min_stability='0.8', existing_output=False, top_n='50', max_features='200', p_max='', q_max='', use_all_samples=False):
    if not counts or not output:
        raise ValueError('请选择丰度表并指定输出目录')
    return dict(counts_path=counts, output_dir=output, metadata_path=metadata or None,
                sample_column=sample_column, group_column=group_column,
                min_prevalence=float(prevalence), min_abundance=float(abundance),
                rho_min=float(rho), permutations=int(permutations), repeats=int(repeats),
                independent=independent, taxonomy_path=taxonomy or None, min_stability=float(min_stability), existing_output=existing_output,
                top_n=int(top_n), max_features=int(max_features),
                p_max=float(p_max) if p_max.strip() else None, q_max=float(q_max) if q_max.strip() else None, use_all_samples=use_all_samples)


def main():
    root = tk.Tk(); root.title(f'小样本微生物统计关联网络探索 v{VERSION}'); root.geometry('900x900')
    frame = ttk.Frame(root, padding=15); frame.pack(fill='both', expand=True)
    values = {k: tk.StringVar(value=v) for k,v in {
        'counts':'', 'metadata':'', 'taxonomy':'', 'output':'', 'sample_column':'sample-id', 'group_column':'group',
        'prevalence':'0.5', 'abundance':'0.001', 'top_n':'50', 'max_features':'200', 'rho':'0.7', 'stability':'0.8', 'p_max':'', 'q_max':'', 'permutations':'9999', 'repeats':'36'}.items()}
    existing_output = tk.BooleanVar(value=False)
    def choose_output():
        directory = filedialog.askdirectory()
        if directory:
            values['output'].set(directory if existing_output.get() else directory + '/network_result')
    for row, (label, key) in enumerate([('丰度表 TSV/CSV','counts'),('元数据（可选）','metadata'),('注释表（可选）','taxonomy'),('输出目录','output'),
                                        ('元数据样本列','sample_column'),('元数据分组列','group_column'),
                                        ('最小出现率（0–1）','prevalence'),('最低平均相对丰度','abundance'),
                                        ('保留丰度前 N 个节点','top_n'),('计算节点上限','max_features'),
                                        ('候选边 |ρ| 阈值','rho'),('最小恢复频率','stability'),
                                        ('原始 p 值小于（可选）','p_max'),('BH q 值小于（可选）','q_max'),
                                        ('随机置换次数（>9 样本）','permutations'),('子样本次数','repeats')]):
        ttk.Label(frame,text=label).grid(row=row,column=0,sticky='w',pady=3)
        ttk.Entry(frame,textvariable=values[key],width=55).grid(row=row,column=1,sticky='ew',padx=8)
        if key in ('counts','metadata','taxonomy'):
            ttk.Button(frame,text='选择',command=lambda k=key: values[k].set(filedialog.askopenfilename() or values[k].get())).grid(row=row,column=2)
        if key == 'output':
            ttk.Button(frame,text='选择目录',command=choose_output).grid(row=row,column=2)
    ttk.Radiobutton(frame,text='新建目录（默认）',variable=existing_output,value=False).grid(row=16,column=0,sticky='w')
    ttk.Radiobutton(frame,text='使用已有目录（拒绝覆盖同名结果）',variable=existing_output,value=True).grid(row=16,column=1,sticky='w')
    ttk.Label(frame,text='新建模式：选择父目录后，请检查并修改末尾的文件夹名称；目标目录必须不存在。',wraplength=820).grid(row=17,column=0,columnspan=3,sticky='w',pady=4)
    use_all_samples = tk.BooleanVar(value=False)
    ttk.Checkbutton(frame,text='合并全部样本构建全局网络（默认按组分别分析）',variable=use_all_samples).grid(row=18,column=0,columnspan=3,sticky='w',pady=4)
    ttk.Label(frame,text=POOLED_NOTICE,wraplength=820,justify='left').grid(row=19,column=0,columnspan=3,sticky='w',pady=3)
    independent = tk.BooleanVar(value=False)
    ttk.Checkbutton(frame,text='我确认当前分析样本独立，且置换所需的可交换性成立',variable=independent).grid(row=20,column=0,columnspan=3,sticky='w',pady=4)
    ttk.Label(frame,text='独立：样本来自独立实验单位；重复采样、配对或技术重复通常不独立。\n可交换：无关联假设下可自由打乱标签；时间、批次等结构可能破坏此条件。本工具不支持受限置换。',wraplength=820,justify='left').grid(row=21,column=0,columnspan=3,sticky='w',pady=4)
    ttk.Label(frame,text=EXPLORATION_NOTICE,wraplength=820,justify='left').grid(row=22,column=0,columnspan=3,sticky='w',pady=4)
    status = tk.StringVar(value='p/q 阈值留空则不筛选；使用 p/q 筛选须确认置换条件。≤9 样本精确置换。')
    ttk.Label(frame,textvariable=status,wraplength=820).grid(row=24,column=0,columnspan=3,sticky='w',pady=8)
    messages = queue.Queue()
    def work(options):
        try: messages.put(('ok', analyze(**options)))
        except Exception as exc: messages.put(('error', str(exc)))
    def poll():
        try: kind, result = messages.get_nowait()
        except queue.Empty:
            root.after(300,poll); return
        button.config(state='normal')
        if kind == 'ok':
            status.set('完成：结果保存在指定目录。'); messagebox.showinfo('分析完成','结果保存在指定的输出目录。')
        else:
            status.set('分析失败：' + result); messagebox.showerror('分析失败',result)
    def start():
        try:
            opt = build_options(*(values[k].get() for k in ('counts','output','metadata','sample_column','group_column','prevalence','abundance','rho','permutations','repeats')),independent.get(),values['taxonomy'].get(),values['stability'].get(),existing_output.get(),values['top_n'].get(),values['max_features'].get(),values['p_max'].get(),values['q_max'].get(),use_all_samples.get())
        except ValueError as exc: messagebox.showerror('输入错误',str(exc)); return
        button.config(state='disabled'); status.set('正在分析，请保持窗口打开……')
        threading.Thread(target=work,args=(opt,),daemon=True).start(); root.after(300,poll)
    button = ttk.Button(frame,text='开始分析',command=start); button.grid(row=23,column=0,sticky='w')
    frame.columnconfigure(1,weight=1); root.mainloop()


if __name__ == '__main__':
    main()
