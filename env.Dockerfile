FROM julia:1.10.0 as julia
FROM python:3.12.1-slim as base

# Julia config
ENV GKSwstype 100
ENV JULIA_CI 'true'
ENV JULIA_NUM_THREADS 'auto'
# Let PythonCall use built-in python
ENV JULIA_CONDAPKG_BACKEND 'Null'
ENV JULIA_PATH '/usr/local/julia/'
ENV JULIA_DEPOT_PATH '/srv/juliapkg/'
ENV PATH ${JULIA_PATH}/bin:${PATH}
COPY --from=julia ${JULIA_PATH} ${JULIA_PATH}

FROM base

WORKDIR /work

# Python dependencies. e.g. matplotlib
RUN pip install --no-cache-dir matplotlib nbconvert

# Julia environment
COPY Project.toml Manifest.toml ./
RUN julia --color=yes -e 'using Pkg; Pkg.add(["IJulia"]); Pkg.activate("."); Pkg.instantiate(); Pkg.precompile()'
