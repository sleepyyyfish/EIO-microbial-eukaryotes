# Microbial co-occurrence network analysis

## Ecological questions addressed by network analysis

Co-occurrence network analysis was performed to investigate three major ecological questions:

### 1. How is the internal structure of microbial communities organized?

Microbial communities can be represented as networks, where:

- **Nodes** represent important ASVs (appear in over 30% samples)
- **Edges** represent significant correlations between ASVs (increase or decrease together,|r|>0.6)

Network topology was evaluated to determine whether microbial communities form weakly connected assemblages or tightly interacting systems.

### 2. Which taxa play key roles in maintaining network organization?

Highly abundant taxa are not necessarily the most ecologically important species. Network analysis identifies key taxa based on their connectivity patterns.

Two centrality measures were calculated:

- **Degree** represents number of connections of each ASV; reflects interaction intensity
High-degree ASVs:
Have many co-occurrence relationships
Participate in multiple microbial interactions
May function as network hubs
- **Betweenness** measures how frequently an ASV appears on the shortest paths connecting other nodes.
High-betweenness ASVs:
Connect different parts of the network
Act as bridge taxa
May influence information or resource flow within the community

### 3. Whether the network is stable or not? 
Network stability was assessed to evaluate how microbial communities respond to species loss or environmental disturbances.

Three complementary approaches were applied:

- 1.**Fragmentation**
Fragmentation evaluates network integrity after removal of important nodes.
Higher fragmentation indicates a more vulnerable network structure.

- 2.**Robustness**
Robustness evaluates resistance to random or targeted species loss.
Two removal strategies were applied: **Random removal**,**Targeted removal**

- 3.**Vulnerability**
Vulnerability was quantified as the maximum decrease in network efficiency caused by the sequential removal of individual nodes, reflecting the dependence of network functioning on specific taxa. 
