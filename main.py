import pandas as pd
from itertools import chain
from itertools import pairwise
from collections import Counter
import pygraphviz as pgv
from IPython.display import Image, display
import ipywidgets as widgets
import os

def format_label_value(val, mode):
  if mode == 'coverage':
    return f"{val:.1f}%"
  else:
    try:
      return str(int(val))
    except Exception:
      return str(val)



df = pd.read_csv('repairExample.csv')
df['start'] = pd.to_datetime(df['Start Timestamp'])
dfs = df[['Case ID', 'Activity', 'start']]

# totals and per-case counts
num_cases = dfs['Case ID'].nunique()

# absolute occurrences (total events)
ev_absolute = dfs['Activity'].value_counts()

# per-case activity counts: dataframe with one row per (case, activity)
case_activity_counts = dfs.groupby(['Case ID', 'Activity']).size().reset_index(name='count')

# case frequency: number of distinct cases containing the activity
ev_case_freq = case_activity_counts.groupby('Activity')['Case ID'].nunique()

# max repetitions: max occurrences of activity within a single case
ev_max_reps = case_activity_counts.groupby('Activity')['count'].max()

# case coverage: percentage of cases containing the activity
ev_case_coverage = (ev_case_freq / num_cases) * 100.0

# helper to get chosen event metric
def get_event_metric_series(mode):
  if mode == 'absolute':
    return ev_absolute
  if mode == 'case':
    return ev_case_freq
  if mode == 'max':
    return ev_max_reps
  if mode == 'coverage':
    return ev_case_coverage

# build per-case traces (tuple of activities per case)
dfs_traces = (dfs.sort_values(by=['Case ID', 'start']).groupby(['Case ID']).agg({'Activity': lambda x: tuple(x)}))

# build per-case pair counters
per_case_pairs = {}
for case_id, row in dfs_traces.iterrows():
  trace = row['Activity']
  pairs = list(pairwise(trace))
  per_case_pairs[case_id] = Counter(pairs)

# aggregate flow metrics
flow_abs = Counter()
flow_case_freq = Counter()
flow_max = Counter()
for case_id, cnt in per_case_pairs.items():
  for pair, c in cnt.items():
    flow_abs[pair] += c
    if c > 0:
      flow_case_freq[pair] += 1
    flow_max[pair] = max(flow_max.get(pair, 0), c)

flow_case_coverage = {pair: (flow_case_freq[pair] / num_cases) * 100.0 for pair in flow_abs.keys()}

def get_flow_metric_dict(mode):
  if mode == 'absolute':
    return dict(flow_abs)
  if mode == 'case':
    return dict(flow_case_freq)
  if mode == 'max':
    return dict(flow_max)
  if mode == 'coverage':
    return dict(flow_case_coverage)

# collect start/end events
ev_start_set = set()
ev_end_set = set()
for case_id, row in dfs_traces.iterrows():
  trace = row['Activity']
  if len(trace) == 0:
    continue
  ev_start_set.add(trace[0])
  ev_end_set.add(trace[-1])
def render_graphs(mode, filter_by='none', thr_ev=0.0, thr_flow=0.0):
  """Render both simple and weighted DFG images using the chosen mode and optional thresholds.

  mode: one of 'absolute','case','max','coverage'
  filter_by: 'none','events','flows','both'
  thr_ev / thr_flow: numeric thresholds (for coverage use percent 0-100)
  """
  event_metric = get_event_metric_series(mode)
  flow_metric = get_flow_metric_dict(mode)

  # build DFG mapping using chosen flow metric
  dfg = {}
  for (a, b), val in flow_metric.items():
    if a not in dfg:
      dfg[a] = Counter()
    dfg[a][b] = val

  # determine color bounds from event_metric
  event_metric_values = [v for v in (event_metric.reindex(index=dfs['Activity'].unique(), fill_value=0).to_dict().values())]
  if len(event_metric_values) > 0:
    color_min = min(event_metric_values)
    color_max = max(event_metric_values)
  else:
    color_min = 0
    color_max = 1

  filter_events = filter_by in ('events', 'both')
  filter_flows = filter_by in ('flows', 'both')

  # Simple event-only visualization
  Gs = pgv.AGraph(strict=False, directed=True)
  Gs.graph_attr['rankdir'] = 'LR'
  Gs.node_attr['shape'] = 'Mrecord'
  all_events = sorted(set(dfg.keys()) | set(d for cs in dfg.values() for d in cs))
  for event in all_events:
    ev_val = event_metric.get(event, 0)
    if filter_events and ev_val < thr_ev:
      continue
    label_val = format_label_value(ev_val, mode)
    Gs.add_node(event, style="rounded,filled", fillcolor="#ffffcc", label=f"{event} ({label_val})")
    if event in dfg:
      for sc in dfg[event]:
        # only add edge if both endpoints pass event threshold (when filtering events)
        if filter_events and event_metric.get(sc, 0) < thr_ev:
          continue
        if filter_flows and dfg[event][sc] < thr_flow:
          continue
        Gs.add_edge(event, sc)
  Gs.draw('simple_heuristic_net.png', prog='dot')
#   display(Image('simple_heuristic_net.png'))

  # Weighted visualization with nodes and edge labels
  trace_counts = sorted(chain(*[c.values() for c in dfg.values()])) if len(dfg) > 0 else [0, 1]
  trace_min = trace_counts[0]
  trace_max = trace_counts[-1]

  Gw = pgv.AGraph(strict=False, directed=True)
  Gw.graph_attr['rankdir'] = 'LR'
  Gw.node_attr['shape'] = 'Mrecord'

  no_traces = num_cases
  # ensure all nodes (including those that only appear as start/end targets) are created
  all_events_full = sorted(set(dfg.keys()) | set(d for cs in dfg.values() for d in cs) | set(ev_start_set) | set(ev_end_set))
  for ev in all_events_full:
    ev_val = event_metric.get(ev, 0)
    # color scale; avoid division by zero
    if color_max != color_min:
      color = int(float(color_min - ev_val) / float(color_min - color_max) * 100.00)
    else:
      color = 5
    my_color = "#ff9933" + str(hex(max(0, min(255, color)))[2:])
    label_val = format_label_value(ev_val, mode)
    Gw.add_node(ev, style="rounded,filled", fillcolor=my_color, label=f"{ev} ({label_val})")

  Gw.add_node("start", shape="circle", label="")
  for ev_start in sorted(ev_start_set):
    if filter_events and event_metric.get(ev_start, 0) < thr_ev:
      continue
    Gw.add_edge("start", ev_start, penwidth=4 * no_traces / (trace_max - trace_min + 1) + 0.1, label=no_traces)

  for event, succesors in dfg.items():
    ev_val = event_metric.get(event, 0)
    if filter_events and ev_val < thr_ev:
      continue
    # color scale; avoid division by zero
    if color_max != color_min:
      color = int(float(color_min - ev_val) / float(color_min - color_max) * 100.00)
    else:
      color = 5
    my_color = "#ff9933" + str(hex(max(0, min(255, color)))[2:])
    label_val = format_label_value(ev_val, mode)
    Gw.add_node(event, style="rounded,filled", fillcolor=my_color, label=f"{event} ({label_val})")
    for succesor, cnt in succesors.items():
      if filter_flows and cnt < thr_flow:
        continue
      if filter_events and event_metric.get(succesor, 0) < thr_ev:
        continue
      lbl = format_label_value(cnt, mode)
      pen = 4 * float(cnt) / (trace_max - trace_min + 1) + 0.1
      Gw.add_edge(event, succesor, penwidth=pen, label=lbl)

  Gw.add_node("end", shape="circle", label="", penwidth='3')
  for ev_end in sorted(ev_end_set):
    if filter_events and event_metric.get(ev_end, 0) < thr_ev:
      continue
    ev_val = event_metric.get(ev_end, 0)
    if color_max != color_min:
      color = int(float(color_min - ev_val) / float(color_min - color_max) * 100.00)
    else:
      color = 0
    my_color = "#ff9933" + str(hex(max(0, min(255, color)))[2:])
    label_val = format_label_value(ev_val, mode)
    Gw.add_node(ev_end, style="rounded,filled", fillcolor=my_color, label=f"{ev_end} ({label_val})")
    Gw.add_edge(ev_end, "end", penwidth=4 * no_traces / (trace_max - trace_min + 1) + 0.1, label=no_traces)

  Gw.draw('simple_heuristic_net_with_events.png', prog='dot')
  display(Image('simple_heuristic_net_with_events.png'))

  # summary
  total_edges = Gw.number_of_edges()
  sum_edge_labels = 0
  for e in Gw.edges():
    try:
      lab = e.attr['label']
      if lab is None:
        continue
      if isinstance(lab, str) and lab.endswith('%'):
        sum_edge_labels += float(lab.rstrip('%'))
      else:
        sum_edge_labels += float(lab)
    except Exception:
      continue
  print(f"Mode: {mode}. Filter: {filter_by}. Thr-events: {thr_ev}. Thr-flows: {thr_flow}. Number of edges: {total_edges}, Label sum: {sum_edge_labels}")

# render once with CLI args (if running as script)
MODE = "absolute"  # Default mode
filter_by = 'none'  # Default filter
threshold_events = 0.0  # Default threshold for events
threshold_flows = 0.0  # Default threshold for flows

# Notebook UI: ipywidgets to select mode and thresholds and regenerate graphs
# Do NOT render graphs on import — show widgets first and render only on button click.
mode_dropdown = widgets.Dropdown(options=['absolute', 'case', 'max', 'coverage'], value=MODE, description='Mode:')
filter_dropdown = widgets.Dropdown(options=['none', 'events', 'flows', 'both'], value=filter_by, description='Filter:')
# dynamically set max values based on data (use .max() for Series)
event_metric = get_event_metric_series(MODE)
try:
  max_ev = int(event_metric.max())
except Exception:
  max_ev = 100 if MODE == 'coverage' else 1
flow_metric = get_flow_metric_dict(MODE)
try:
  max_flow = int(max(flow_metric.values())) if len(flow_metric) > 0 else 1
except Exception:
  max_flow = 1
# coverage mode expects percent thresholds
if MODE == 'coverage':
  max_ev = 100
  max_flow = 100

thr_ev_box = widgets.IntSlider(value=min(int(threshold_events), max_ev), min=0, max=max_ev if max_ev>0 else 1, step=1, description='Thr events:')
thr_flow_box = widgets.IntSlider(value=min(int(threshold_flows), max_flow), min=0, max=max_flow if max_flow>0 else 1, step=1, description='Thr flows:')

# update slider maxima when mode dropdown changes
def _on_mode_change(change):
  if change.get('name') != 'value':
    return
  new_mode = change.get('new')
  ev_m = get_event_metric_series(new_mode)
  try:
    new_max_ev = int(ev_m.max())
  except Exception:
    new_max_ev = 100 if new_mode == 'coverage' else 1
  flow_m = get_flow_metric_dict(new_mode)
  try:
    new_max_flow = int(max(flow_m.values())) if len(flow_m) > 0 else 1
  except Exception:
    new_max_flow = 1
  if new_mode == 'coverage':
    new_max_ev = 100
    new_max_flow = 100
  # set slider maxima and clamp values if needed
  thr_ev_box.max = new_max_ev if new_max_ev > 0 else 1
  thr_flow_box.max = new_max_flow if new_max_flow > 0 else 1
  if thr_ev_box.value > thr_ev_box.max:
    thr_ev_box.value = thr_ev_box.max
  if thr_flow_box.value > thr_flow_box.max:
    thr_flow_box.value = thr_flow_box.max

mode_dropdown.observe(_on_mode_change, names='value')
generate_btn = widgets.Button(description='Generate graph')
out = widgets.Output()

def _on_generate_clicked(b):
  # clear previous output (graphs) but keep widgets visible
  out.clear_output(wait=True)
  mode_w = mode_dropdown.value
  filter_w = filter_dropdown.value
  thr_ev_w = float(thr_ev_box.value)
  thr_flow_w = float(thr_flow_box.value)
  # remove old graph files if present to avoid stale images
  for fn in ('simple_heuristic_net.png', 'simple_heuristic_net_with_events.png'):
    try:
      if os.path.exists(fn):
        os.remove(fn)
    except Exception:
      pass
  with out:
    render_graphs(mode_w, filter_w, thr_ev_w, thr_flow_w)

generate_btn.on_click(_on_generate_clicked)

ui = widgets.HBox([mode_dropdown, filter_dropdown, thr_ev_box, thr_flow_box, generate_btn])
display(ui)
display(out)