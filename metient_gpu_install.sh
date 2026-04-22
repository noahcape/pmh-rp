# 1. Create conda/mamba environment
mamba create -n metient_gpu_2 -c conda-forge python=3.9 -y
mamba activate metient_gpu_2

# 2. Install PyTorch 1.13 with CUDA 11.7 (supports P100)
mamba install pytorch=1.13.1 torchvision torchaudio cudatoolkit=11.7 -c pytorch -y

# 3. Install other dependencies from conda-forge
mamba install -c conda-forge \
    numpy=1.24 \
    scipy=1.10 \
    matplotlib=3.7 \
    pandas=1.5 \
    seaborn=0.12 \
    joblib=1.3 \
    networkx=3.1 \
    ipython \
    ipykernel \
    graphviz \
    pygraphviz \
    pillow=10.1 \
    tqdm \
    pydot \
    -y

# 4. Install metient itself via pip
pip install --no-deps metient
