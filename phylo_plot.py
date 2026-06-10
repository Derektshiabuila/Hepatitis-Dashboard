from Bio import Phylo
import io
import plotly.graph_objects as go

def build_tree_figure(newick_str: str) -> go.Figure:
    """Parse Newick string and return a Plotly rectangular tree figure."""
    if not newick_str:
        return go.Figure()

    try:
        tree = Phylo.read(io.StringIO(newick_str), "newick")
    except Exception as e:
        fig = go.Figure()
        fig.add_annotation(text=f"Error parsing tree: {e}", showarrow=False)
        return fig
        
    tree.ladderize()
    
    y_coord = {}
    y_current = 0
    
    # Assign y coordinates (terminals first, then internal nodes average)
    for cl in tree.find_clades(order='postorder'):
        if cl.is_terminal():
            y_coord[cl] = y_current
            y_current += 1
        else:
            if cl.clades:
                y_coord[cl] = sum(y_coord[c] for c in cl.clades) / len(cl.clades)
            else:
                y_coord[cl] = y_current
            
    x_coord = {}
    x_coord[tree.root] = 0
    for cl in tree.find_clades(order='preorder'):
        for child in cl.clades:
            bl = child.branch_length if child.branch_length else 0
            x_coord[child] = x_coord[cl] + bl

    lines_x = []
    lines_y = []
    
    text_x = []
    text_y = []
    text_labels = []
    text_colors = []
    
    for cl in tree.find_clades(order='preorder'):
        if cl.clades:
            y_children = [y_coord[c] for c in cl.clades]
            if y_children:
                ymin, ymax = min(y_children), max(y_children)
                # Vertical line
                lines_x.extend([x_coord[cl], x_coord[cl], None])
                lines_y.extend([ymin, ymax, None])
                
                # Horizontal lines to children
                for child in cl.clades:
                    lines_x.extend([x_coord[cl], x_coord[child], None])
                    lines_y.extend([y_coord[child], y_coord[child], None])
                
        if cl.name:
            text_x.append(x_coord[cl])
            text_y.append(y_coord[cl])
            text_labels.append(cl.name)
            is_user = not cl.name.startswith("ref_")
            text_colors.append("#d9534f" if is_user else "#6c757d")

    fig = go.Figure()
    
    # Add lines
    fig.add_trace(go.Scatter(
        x=lines_x, y=lines_y,
        mode='lines',
        line=dict(color='#adb5bd', width=1.5),
        hoverinfo='skip'
    ))
    
    # Add nodes
    fig.add_trace(go.Scatter(
        x=text_x, y=text_y,
        mode='markers+text',
        marker=dict(size=8, color=text_colors, line=dict(width=1, color='white')),
        text=text_labels,
        textposition="middle right",
        textfont=dict(size=11, color=text_colors),
        hoverinfo='text'
    ))
    
    # Estimate height based on number of terminals (approx 20px per terminal)
    height = max(400, y_current * 20)
    
    fig.update_layout(
        height=height,
        showlegend=False,
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin=dict(l=10, r=120, t=10, b=10),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, autorange='reversed'),
    )
    
    return fig
