# GPUMDkit Agent Skill

GPUMDkit Agent Skill provides instructions for AI agents working with GPUMDkit
tools, GPUMD simulations, and NEP training and analysis. Once it is installed,
describe what you need: the agent can consult the guidance, choose suitable
commands, and help carry out the task.

Before using it, make sure GPUMDkit is installed and `gpumdkit.sh` runs in the
environment where the agent works. The agent must also be allowed to run shell
commands there. Choose the English or Chinese skill according to your language;
both provide the current GPUMDkit guidance.

You can copy this prompt to your agent:

```text
Run `gpumdkit.sh -skill`. Follow the guidance it prints to install the relevant
GPUMDkit skills in this client's user-level global skills directory, then confirm
that the skills are available.
```

`gpumdkit.sh -skill` displays the skill locations and installation instructions;
the agent performs the installation.

After installation, describe your task directly in natural language. For example:

- “Convert this VASP structure to extxyz and check the result.”
- “Inspect my GPUMD outputs and plot the diffusion results.”

If scientific parameters are unclear, the agent should ask before choosing them;
long-running or expensive calculations require your confirmation first.

## Example Requests

Ask in plain language and name the files. The agent routes the request to the right GPUMDkit commands and shows them before running anything that changes files or takes time. Requests that work well:

- “List the compositions in `train.xyz`, report the energy and force ranges, and tell me if any two atoms are closer than 1 Å.” (analysis only; nothing is modified)
- “Convert the OUTCAR files under `scf/` to extxyz and check the result before I train.”
- “Plot the MSD from `msd.out` and tell me whether the slope looks converged.”
- “Explain what the `ensemble` line in my `run.in` does.”

The agent should stop and ask you in these cases: a scientific choice is open (a temperature, a cutoff, a distance threshold, which species is mobile), a command would overwrite or delete files, or a step would start a simulation, a training or a scheduler job. Answer with the value you want, or tell it to leave the choice to you.
