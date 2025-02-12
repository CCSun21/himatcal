# Development Environment Instructions

Please follow these guidelines for managing development environments:

## Requirements

- micromamba for environment management
- uv for fast Python package installation

## Environment Setup

Use micromamba to create and activate the environment:

```bash
micromamba create -n himatcal
micromamba activate himatcal
```

## Package Management

Use uv for dependency management:

```bash
uv pip install <package-name>
```

For development installations:

```bash
uv pip install -e .
```

## Project Structure

- Keep dependencies minimal and well-documented
- Follow Python best practices for package structure

## Best Practices

1. Always work within the micromamba environment
2. Document new dependencies in `pyproject.toml`
3. Leverage uv's speed for package management