# Microbial network visualization using Gephi
## Workflow summary

```
nodes.csv + edges.csv

        ↓

Import spreadsheet

        ↓

Undirected network

        ↓

Taxonomic coloring:
    - Phylum

        ↓

Edge visualization:
    - Positive correlation: red
    - Negative correlation: gray/black

        ↓
Node visualization:
    - Abundance
    - Degree
    - Betweenness

        ↓
ForceAtlas2 layout

        ↓

Export final network figure
```

## 1. Input files

The network visualization was performed using **Gephi**.

Two files were imported:

- `nodes.csv`
    - Node information table
    - Contains ASV IDs and node attributes, such as phylum taxonomy and abundance

Example:

|Id|Phylum|Abundance|
|ASV1|Chlorophyta|133726|
|ASV2|Cercozoa|65969|


- `edges.csv`
    - Network edge table
    - Contains interactions between ASVs
    - Includes edge weight representing correlation strength
    - Positive and negative correlations were retained

Example:

|source|target|weight|pvalue|
|ASV1|ASV2|133726|-0.606628184|0.0076|

---

# 2. Import network files

## Import nodes table

Path:

```
File → Import Spreadsheet
```

Settings:

- Import as: **Nodes table**
- Graph type: **Undirected**
- Ensure the ID column is recognized as node ID


## Import edges table

Path:

```
File → Import Spreadsheet
```

Settings:

- Import as: **Edges table**
- Graph type: **Undirected**
- Weight column: correlation value


After importing both tables, check that the network structure is correctly displayed.

---

# 3. Node color by taxonomy

Node colors were assigned according to phylum classification.

Path:

```
Appearance → Nodes → Partition → Color
```

Select attribute:

```
Phylum
```

Purpose:

- Display taxonomic composition of the microbial network
- Compare distribution of different microbial groups


---

# 4. Edge color visualization

Edge color was adjusted according to correlation direction.

Path:

```
Appearance → Edges → Weight
```


Settings:
(black--red)
- Negative correlations:
    - gray/black edges

- Positive correlations:
    - red edges


Edge thickness keeps the same


---

# 5. Node size visualization

Node size was adjusted according to different network properties.

Path:

```
Appearance → Nodes → Ranking → Size
```


Three types of node ranking were generated:


## Node size by abundance

Select:

```
Ranking → Abundance → Size
```

Purpose:

- Display dominant ASVs with larger nodes
- Reflect relative importance based on abundance


---

## Node size by degree

Select:

```
Ranking → Degree → Size
```

Purpose:

- Highlight highly connected ASVs
- Identify potential network hubs


---

## Node size by betweenness centrality

Before visualization:

```
Statistics → Network Diameter
```

Calculate network statistics to obtain betweenness values.

Then:

```
Appearance → Nodes → Ranking → Betweenness → Size
```

Purpose:

- Highlight ASVs acting as bridges between different network modules
- Identify key connector taxa


---

# 6. Network layout

The network was visualized using:

```
Layout → ForceAtlas2
```

Parameters can be adjusted according to network complexity.

ForceAtlas2 was used to generate a force-directed layout, where highly connected nodes tend to cluster together and network structures become visually interpretable.


---

# 7. Export network figure

After adjusting:

- Layout
- Node size
- Node color
- Edge color
- Background


Preview:

```
Preview → Refresh
```


Export:

```
File → Export
```

Recommended format:

- SVG
- PDF
- PNG
