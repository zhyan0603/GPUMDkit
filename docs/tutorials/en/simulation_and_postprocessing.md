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
