# Archipelago Explorer

**Frequency, Persistence and Regime Visualization Toolkit**

Archipelago Explorer is an exploratory visualization tool. It does not claim
to discover new particles, new laws, or hidden physical symmetries. It provides
an interface for testing whether frequency, lifetime, persistence cycles, and
phase mappings reveal robust structures in known physical datasets.

## Run With Docker

```powershell
docker compose up --build -d
```

Open [http://localhost:8000](http://localhost:8000).

Stop the app with:

```powershell
docker compose down
```

## Local Development

```powershell
npm install
npm run dev
npm test
npm run build
```

## Features

- CSV upload and dataset preview
- Stable-particle handling controls
- Standard and exploratory formula catalog
- Preset and custom logarithmic bases
- 2D and 3D Plotly visualizations
- Base-sweep animation
- Gap and circular metrics
- Persistent-neighbor network
- Shuffled and randomized controls
- CSV, JSON, Markdown, PNG, and SVG exports

## License

License: MIT. See [../LICENSE](../LICENSE).

Read [the scientific scope](docs/scientific_scope.md) before interpreting
visual structures.
