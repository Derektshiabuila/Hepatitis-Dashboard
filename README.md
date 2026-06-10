# Hepatitis Virus Sequence Dashboard

A web-based sequence analysis and epidemiological visualization dashboard for Hepatitis B (HBV), Hepatitis C (HCV), and Hepatitis E (HEV) viruses. 

The dashboard consists of two major parts:
1. **Interactive Global Analytics**: Visualization of WHO regional statistics, genotype distributions, temporal trends, and country-level burden data (sourced from IHME GBD and WHO).
2. **User Sequence Analysis ("My Sequences")**: Fast end-to-end bioinformatics pipeline to:
   - Identify genotypes and nearest reference sequences (via MAFFT and EPA-ng phylogenetic placement).
   - Identify drug resistance profiles and substitutions of interest (via GLUE).
   - Reconstruct phylogenetic subtrees (via IQ-TREE2).
   - Detect recombinant events (via RDP5 running under Wine).

---

## Getting Started: Local Installation

### Prerequisites
1. **Conda**: Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda.
2. **Docker**: Required if you run GLUE mutation analysis (which runs inside containers).
3. **Wine**: Required to run RDP5CL.exe (on macOS/Linux).

### Setup
1. Clone this repository:
   ```bash
   git clone https://github.com/Derektshiabuila/Hepatitis-Dashboard.git
   cd Hepatitis-Dashboard
   ```
2. Provision the Conda environment:
   ```bash
   conda env create -f envs/phylo.yaml
   conda activate phylo
   ```
3. Install Python dependencies:
   ```bash
   pip install -r requirements.txt
   ```
4. Start the dashboard:
   ```bash
   python Full_Hepatitis_page.py
   ```
   Open `http://127.0.0.1:8051` in your browser.

---

## Running with Docker

To make hosting the dashboard simple, a production-ready `Dockerfile` is provided. It sets up Debian, installs Wine, Conda, Python dependencies, and provisions the `phylo` tools environment.

### 1. Build the Docker Image
```bash
docker build -t hepatitis-dashboard .
```

### 2. Run the Container
Because the sequence analysis pipeline launches GLUE inside docker containers, you must mount the host's Docker socket into the dashboard container (**Docker-out-of-Docker** mode). This allows the dashboard to launch sibling containers on the host.

```bash
docker run -d \
  -p 8051:8051 \
  -v /var/run/docker.sock:/var/run/docker.sock \
  --name hep-dash \
  hepatitis-dashboard
```

---

## Cloud Deployment (e.g. Render)

To host the dashboard publicly so anyone can access it:

1. **Create an account** on [Render.com](https://render.com) using your GitHub profile.
2. Click **New +** > **Web Service**.
3. Connect the `Hepatitis-Dashboard` repository.
4. Select the **Docker** runtime (instead of Python). Render will automatically detect the `Dockerfile` in the root of your repository.
5. Select the **Free** tier (or paid tiers if you require persistent storage or higher resources for heavy alignment tasks).
6. Click **Deploy Web Service**. Within a few minutes, Render will provision the container and provide you with a live URL (e.g. `https://hepatitis-dashboard.onrender.com`).

*Note: Headless Wine and Conda tools (MAFFT/EPA-ng/RDP5) will work fully on cloud environments. If you want GLUE database containers running in the cloud, you should deploy to a VPS (e.g., AWS EC2 or DigitalOcean) where you can run a Docker Compose stack.*
