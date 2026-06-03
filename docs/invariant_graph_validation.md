# Invariant Graph Validation

Persistent graph structure does not imply physical interaction. It indicates only that selected entities remain close under the chosen mathematical representations and tolerances.

All outputs are exploratory unless they pass randomized controls and tolerance checks. No new physics or discovery is claimed.

## Validation Settings
- Particles/entities tested: 29
- Bases tested: 2009
- Randomized controls use the generated representative control-base grid recorded in `graph_validation_summary.csv`.
- Monte Carlo runs per control: 1000
- Primary angular tolerance: 10 degrees
- Primary persistence threshold: 0.7
- Empirical p-value threshold: 0.05

## Conclusion So Far

The current validation run finds a small set of persistent graph edges under the primary angular tolerance of 10 degrees across 2009 tested bases. The clearest base-persistent angular examples are B0 / Bs0, W / Z, B_plus / Bs0, and omega_782 / top, each reported with explicit persistence fractions and empirical p-values in `data/invariant_edges.csv`.

This is not yet evidence for a physical relation. Several robust-looking clusters in distance and k-nearest-neighbor graphs have high persistence but do not pass randomized-control p-value checks, so they should be treated as weak candidates or likely artifacts until stronger controls support them.

Family modularity is the most notable label-level signal in this run for the angular graph, but it remains tolerance- and control-dependent. Interaction labels are less clearly predictive in the current summary. The observed structure is still strongly tied to the selected N-derived phase representation, so the safe interpretation is exploratory mathematical persistence, not physical interaction.

The practical conclusion is: use the persistent edges and clusters as candidates for further robustness testing, not as discoveries. Results that disappear under N shuffling, random phase controls, or tolerance changes should be labeled tolerance-sensitive, representation-sensitive, or likely artifacts.

## 1. Which relationships survive changes of base?
- B0 / Bs0: tolerance angular_close_deg=10.0, 1985/2009 bases, persistence 0.988, empirical p=0.0264.
- B0 / Bs0: tolerance angular_close_deg=10.0, 1985/2009 bases, persistence 0.988, empirical p=0.0244.
- W / Z: tolerance angular_close_deg=10.0, 1771/2009 bases, persistence 0.882, empirical p=0.0222.
- W / Z: tolerance angular_close_deg=10.0, 1771/2009 bases, persistence 0.882, empirical p=0.0200.
- B_plus / Bs0: tolerance angular_close_deg=10.0, 1684/2009 bases, persistence 0.838, empirical p=0.0266.
- B_plus / Bs0: tolerance angular_close_deg=10.0, 1684/2009 bases, persistence 0.838, empirical p=0.0294.
- omega_782 / top: tolerance angular_close_deg=10.0, 1654/2009 bases, persistence 0.823, empirical p=0.0182.
- omega_782 / top: tolerance angular_close_deg=10.0, 1654/2009 bases, persistence 0.823, empirical p=0.0390.

## 2. Which relationships survive changes of representation?
The multi-representation graph requires closeness in at least three representations. Relationships absent from that graph are labeled representation-sensitive rather than robust.

## 3. Which relationships survive randomized controls?
70 graph metric rows pass p < 0.05. Rows that do not pass should be treated as weak candidates or likely artifacts.

## 4. Which relationships are tolerance-sensitive?
Mean persistence by reported tolerance: {"2.0": 1.0, "10.0": 0.1332}. Strong changes across this map indicate tolerance-sensitive structure.

## 5. Which graph metrics are stronger than random?
- angular_closeness, average_shortest_path, tolerance angular_close_deg=10.0: observed 1.2, null mean 0.822, null std 0.4038, empirical p=0.0440.
- angular_closeness, modularity_by_family, tolerance angular_close_deg=10.0: observed 0.2812, null mean -0.1526, null std 0.2386, empirical p=0.0250.
- angular_closeness, average_shortest_path, tolerance angular_close_deg=10.0: observed 1.2, null mean 0.8389, null std 0.3829, empirical p=0.0330.
- angular_closeness, modularity_by_family, tolerance angular_close_deg=10.0: observed 0.2812, null mean -0.1149, null std 0.2345, empirical p=0.0400.
- angular_closeness, average_edge_weight, tolerance angular_close_deg=10.0: observed 0.8828, null mean 0.88, null std 7.66e-17, empirical p=0.0010.
- angular_closeness, modularity_by_family, tolerance angular_close_deg=10.0: observed 0.2812, null mean -0.09372, null std 0.1683, empirical p=0.0410.
- angular_closeness, number_of_edges, tolerance angular_close_deg=10.0: observed 4, null mean 1.035, null std 0.9674, empirical p=0.0140.
- angular_closeness, average_shortest_path, tolerance angular_close_deg=10.0: observed 1.2, null mean 0.671, null std 0.478, empirical p=0.0210.
- angular_closeness, degree_distribution, tolerance angular_close_deg=10.0: observed 0.5184, null mean 0.2033, null std 0.157, empirical p=0.0170.
- angular_closeness, largest_component_size, tolerance angular_close_deg=10.0: observed 3, null mean 1.699, null std 0.5276, empirical p=0.0350.
- angular_closeness, modularity_by_family, tolerance angular_close_deg=10.0: observed 0.2812, null mean -0.1504, null std 0.2318, empirical p=0.0160.
- angular_closeness, average_shortest_path, tolerance angular_close_deg=10.0: observed 1.2, null mean 0.8536, null std 0.3738, empirical p=0.0420.

## 6. Which clusters are most robust?
- Cluster 2: tau;D0;D_plus;Ds_plus;B_plus;B0;Bs0; tolerance euclidean_close_percentile=10.0; persistence 1.000; family purity 0.857; interaction purity 1.000; p=0.5111.
- Cluster 9: J_psi;Upsilon_1S; tolerance euclidean_close_percentile=10.0; persistence 1.000; family purity 1.000; interaction purity 1.000; p=0.5019.
- Cluster 7: pi_plus;K_plus;K_long; tolerance euclidean_close_percentile=10.0; persistence 1.000; family purity 1.000; interaction purity 1.000; p=0.5076.
- Cluster 4: Lambda;Sigma_plus;Xi_minus;Omega_minus;K_short; tolerance euclidean_close_percentile=10.0; persistence 1.000; family purity 0.800; interaction purity 1.000; p=0.5176.
- Cluster 10: J_psi;Upsilon_1S; tolerance angular_close_deg=10.0; persistence 1.000; family purity 1.000; interaction purity 1.000; p=0.3163.
- Cluster 6: B_plus;B0;Bs0; tolerance neighbor_rank_k=2.0; persistence 1.000; family purity 1.000; interaction purity 1.000; p=0.4560.
- Cluster 7: W;Z;top; tolerance neighbor_rank_k=2.0; persistence 1.000; family purity 0.667; interaction purity 1.000; p=0.5630.
- Cluster 10: W;Z;top; tolerance euclidean_close_percentile=10.0; persistence 1.000; family purity 0.667; interaction purity 1.000; p=0.5119.

## 7. Are family or interaction labels predictive?
Use `modularity_by_family` and `modularity_by_interaction` in `data/graph_validation_summary.csv`. A label is treated as predictive only when observed modularity exceeds the matching randomized-label controls with p < 0.05.

## 8. Is the observed structure mainly driven by N?
Controls that shuffle or randomize N directly test this. If structure disappears under `shuffle_N_values`, `random_logN_same_range`, or `random_phase_same_z`, the reported graph is mainly driven by the chosen N projection.

## 9. Is there evidence beyond visualization?
Only rows with explicit tolerance, base count, Monte Carlo count, observed metric, null mean, null std, and empirical p-value provide evidence beyond visualization. Non-significant rows are reported as non-significant.

## 10. What should not be claimed?
Do not claim discovery, new interaction, preferred mathematical base, or physical causation. Persistent closeness is a property of selected representations, tolerances, and data preprocessing.
