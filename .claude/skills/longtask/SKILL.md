---
name: longtask
description: Hand a long-running shell command to the user as a self-contained backgrounded job — taskset core-pinned, BLAS-thread capped, redirected to a log file in logs/, and ALSO persisted as a runnable txt file under txt_task/<descriptive>_cmds.txt. INVOKE ONLY WHEN THE USER EXPLICITLY ASKS via the slash form /longtask, or says "use longtask" / "use the longtask skill". Do NOT invoke proactively for ordinary long commands — Claude's own background-task tools (Bash run_in_background) handle those.
---

# longtask — backgrounded shell command + persisted txt task file

Use this skill ONLY when the user explicitly invokes `/longtask` or asks for "the longtask pattern" / "use longtask". Do not auto-invoke for any long task.

## Two outputs every time

Every invocation must produce BOTH of:

1. A fenced bash block in the chat that the user can copy-paste.
2. A text file at `txt_task/<descriptive_name>_cmds.txt` with the same command plus pre-flight CPU check, monitoring kit, runtime estimate, and follow-up commands. Use the Write tool to create it.

## Pre-flight: pick a free core (filtered to user's allowed set)

Always include a CPU-availability check above the run command (in both the chat reply and the txt file). The machine is shared and the user's account is cgroup-restricted to a specific subset of CPUs (e.g. `112-119,240-247`). Pinning to a core outside that set fails; pinning to a busy core wastes time.

The user's allowed cores come from `taskset -pc $$`. Filter the per-core usage to those. Note: `mpstat` is NOT installed on this machine; use `top -bn2`.

```bash
# 1. Show YOUR allowed cores (from cgroup / affinity)
ALLOWED=$(taskset -pc $$ | awk -F'list: ' '{print $2}')
echo "your allowed cores: $ALLOWED"

# 2. Per-allowed-core usage from top, busiest first.
#    Cores NOT in the output are 100% idle (top suppresses stable rows).
top -bn2 -d 1 -1 2>/dev/null | awk -v allowed="$ALLOWED" '
  BEGIN {
    n=split(allowed, parts, ",")
    for (i=1;i<=n;i++) {
      if (split(parts[i], r, "-") == 2) for (c=r[1]+0;c<=r[2]+0;c++) ok[c]=1
      else ok[parts[i]+0]=1
    }
    snap=0
  }
  /^top - / { snap++ }
  snap==2 && /^%Cpu[0-9]/ {
    c=$1; gsub(/^%Cpu/,"",c); gsub(/[^0-9]/,"",c)
    if ((c+0) in ok) {
      for (i=1;i<=NF;i++) if ($i=="id,") { idle=$(i-1); break }
      printf "cpu %3d  used=%5.1f%%  idle=%5.1f%%\n", c, 100-idle, idle
    }
  }
' | sort -k3 -t= -nr
```

Then plug a low-usage (or absent — those are 100% idle) core into `taskset -c <CORE>` below.

## Command shape

```bash
mkdir -p logs && env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
nohup taskset -c <CORE> ./build/<binary> \
  <flags> \
  > logs/<descriptive_name>.log 2>&1 &
```

Required pieces:

- `mkdir -p logs && ...` — guarantee the log directory exists in one shot.
- `env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1` — pin BLAS to single-threaded so it does not fight `taskset`.
- `nohup` — set SIGHUP-ignore at exec time so the job survives SSH disconnect / login-shell exit. The user runs over SSH from a Mac that may sleep or close; without `nohup`, the shell's SIGHUP on exit kills the job.
- `taskset -c <CORE>` — default core `112` unless the user has been using a different one (check `txt_task/step4_per_step_checklist_cmds.txt` or recent transcripts for the convention; pick a free core if running multiple tasks in parallel).
- `> logs/<name>.log 2>&1 &` — redirect both streams and background. The redirect also prevents `nohup` from creating a `nohup.out` file.

## Naming convention

Same name for the log file and the txt task file (without the `_cmds` suffix on the log).

Encode the experiment identity:
- binary kind (`realtime`, `static`, `imag_time`, …)
- particle number `N{n}` and basis size `K{k}` if relevant
- distinguishing parameter (e.g. `T2pi`, `dt1e3`, `rcond1e-6`, `no_wiener`)

Examples:
- `txt_task/realtime_N2_K5_cmds.txt`         ↔ `logs/realtime_N2_K5.log`
- `txt_task/static_N1_K5_cmds.txt`           ↔ `logs/static_N1_K5.log`
- `txt_task/realtime_N1_K5_no_wiener_cmds.txt` ↔ `logs/realtime_N1_K5_no_wiener.log`

## Concurrent runs that share fixed output paths: sandbox the cwd

The ECG verify binary writes its CSVs through `out_path()` (see `src/verify/common.cpp`), which is hard-coded to the **relative** directory `"out/verify"`. There is **no CLI flag** to redirect it. So if two jobs with the same `(N, K)` run from the same cwd, their final CSVs collide and the second writer wins.

When the user wants to run a variant (different `--rcond`, `--no-wiener`, different solver knob) **concurrently** with an already-running job at the same `(N, K)`, do not try to coordinate output names. Instead, **sandbox the cwd** so each job has its own `out/verify/` namespace.

### Recipe

1. **Build the sandbox cwd** under `runs/<variant_tag>/out/verify/`. Tag examples: `no_wiener_N1K5`, `rcond1e-7_N2K7`.

2. **Symlink any input caches the run needs** into the sandbox's `out/verify/`. For Step-4 that means the Step-1 basis cache (`step1_ecg_basis_N{N}_K{K}.csv`). The binary reads it via the same cwd-relative `out_path()` and a symlink is followed transparently.

3. **Use absolute paths** for the binary and the log redirect, so the `cd` into the sandbox does not break either:

```bash
cd /home/gyqyan/zexuan/ecg1d_cpp

mkdir -p runs/<variant_tag>/out/verify logs
ln -sf /home/gyqyan/zexuan/ecg1d_cpp/out/verify/<input_cache>.csv \
       runs/<variant_tag>/out/verify/<input_cache>.csv

env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
nohup taskset -c <CORE> bash -c '
  cd /home/gyqyan/zexuan/ecg1d_cpp/runs/<variant_tag>
  /home/gyqyan/zexuan/ecg1d_cpp/build/<binary> \
    <flags> \
    <variant flag, e.g. --no-wiener>
' > logs/<name>.log 2>&1 &
echo "PID=$!"
```

The `bash -c` wrapper is required so the `cd` runs **inside** the backgrounded process; the outer shell's cwd does not change. `nohup` is on the outer command — it propagates to the inner binary because the wrapper is not a login shell.

### What lands where

| Job (started from)                    | Output dir                                  |
|---------------------------------------|---------------------------------------------|
| project root                          | `out/verify/step4_ecg_*_N{N}_K{K}.csv`      |
| `runs/<variant_tag>/`                 | `runs/<variant_tag>/out/verify/step4_ecg_*_N{N}_K{K}.csv` |

Logs always land under the project's `logs/` because the redirect path is absolute.

### Plotting from a sandbox result

The plotters (`scripts/plot_step4_*.py`) read the default `out/verify/` paths. To plot a sandbox variant, swap its CSVs into the default slot **only after backing up whatever is there** (typically into `out/verify/<other_variant>/`), plot, then restore. The variant tag in the sandbox path doubles as the archive subdirectory name. Document this swap in the txt file's "When the task finishes" section so the user does not have to figure it out from scratch later.

### When NOT to sandbox

If only one job at a given `(N, K)` is running and you do not need to preserve a previous result, the plain command shape is fine — overwriting `out/verify/*` is the expected behavior.

## Txt file structure

Use this template for `txt_task/<name>_cmds.txt`:

```
Task: <one-line description>

Purpose
- bullet list explaining what the run is for and why these parameters

Pre-flight: pick a free core (filtered to your allowed set)

# 1. Your allowed cores (cgroup / affinity)
ALLOWED=$(taskset -pc $$ | awk -F'list: ' '{print $2}')
echo "your allowed cores: $ALLOWED"

# 2. Per-allowed-core usage; cores not listed are 100% idle.
top -bn2 -d 1 -1 2>/dev/null | awk -v allowed="$ALLOWED" '
  BEGIN {
    n=split(allowed, parts, ",")
    for (i=1;i<=n;i++) {
      if (split(parts[i], r, "-") == 2) for (c=r[1]+0;c<=r[2]+0;c++) ok[c]=1
      else ok[parts[i]+0]=1
    }
    snap=0
  }
  /^top - / { snap++ }
  snap==2 && /^%Cpu[0-9]/ {
    c=$1; gsub(/^%Cpu/,"",c); gsub(/[^0-9]/,"",c)
    if ((c+0) in ok) {
      for (i=1;i<=NF;i++) if ($i=="id,") { idle=$(i-1); break }
      printf "cpu %3d  used=%5.1f%%  idle=%5.1f%%\n", c, 100-idle, idle
    }
  }
' | sort -k3 -t= -nr

How to run

cd /home/gyqyan/zexuan/ecg1d_cpp

mkdir -p logs && env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
nohup taskset -c <CORE> ./build/<binary> \
  <flags> \
  > logs/<name>.log 2>&1 &

Monitor

# live progress
tail -f logs/<name>.log

# is it still alive?
ps -o pid,etime,pcpu,cmd -p $(pgrep -f "<binary> <distinctive flag>") 2>/dev/null

# any blowups or warnings?
grep -E "BLOWUP|WARN|FATAL|drift" logs/<name>.log

Estimated runtime
- give a wall-clock range based on prior runs and the dominant cost scaling

Outputs
- list every artifact the run produces with its path

When the task finishes, follow up with

<plot / analysis commands the user will want next>
```

## Rescue command for an already-running job

If the user already started a long task WITHOUT `nohup` and now wants to close their SSH, give them this one-liner to run in the SAME shell that launched it:

```bash
disown -h %1            # or by PID: disown -h <PID>
```

`disown -h` marks the running job to ignore SIGHUP on shell exit — equivalent to having prefixed the command with `nohup`. The job keeps running and the log keeps being written.

## After-running checklist for the chat reply

In the chat (separate from the txt file), include the same monitoring kit so the user does not need to open the txt file just to start watching:

```bash
tail -f logs/<name>.log
ps -o pid,etime,pcpu,cmd -p $(pgrep -f "<binary> <distinctive flag>") 2>/dev/null
```

And mention the txt file path so the user knows where the durable copy lives.

## Estimate runtime

Always give a rough wall-clock estimate, based on:
- prior runs of the same binary at smaller parameters (cite if available)
- dominant cost scaling — for ECG+TDVP that is roughly `K^2 · N!` per RK4 stage

Range form is fine ("25–40 minutes") when the estimate is uncertain.

## After the task finishes

When the user reports the task is done — or when you check the log and see success markers — proceed with the next thing they asked for, usually plotting or analysis. Do not re-run the simulation in the foreground.

## Anti-patterns

- Do not use Claude's own `run_in_background: true` Bash flag for this skill. The user owns the process.
- Do not run a long task in the foreground "to test it works" before handing the command over.
- Do not omit the log redirect; the user should never have to scroll terminal scrollback to find output.
- Do not skip the env BLAS-thread caps — without them, `taskset` to one core lets BLAS saturate other cores anyway.
- Do not skip writing the txt file — it is the durable record the user comes back to weeks later.
- Do not omit `nohup`. The user runs over SSH and may close their laptop — without `nohup`, the SIGHUP on shell exit kills the job. (User explicitly does not want a tmux/screen wrapper, so `nohup` is the chosen mechanism.)
- Do not use `&>` or `>&` as the redirect — use the portable `> logs/X.log 2>&1` form.
- Do not run two jobs that share the same fixed `(N, K)` output paths concurrently from the same cwd. The binary has no `--out` flag — collision is silent and the last writer wins. Sandbox the cwd per the section above.
- Do not try to "fix" the collision by adding a `--csv-suffix` flag to the binary or by post-renaming inside a single cwd race. The cwd-sandbox is the simple, correct mechanism — both jobs run unchanged and the kernel keeps their writes disjoint.

## Worked example

User: `/longtask realtime N=4 K=7 T=2π`

Claude:

> Wrote `txt_task/realtime_N4_K7_cmds.txt`. Run:
>
> ```bash
> cd /home/gyqyan/zexuan/ecg1d_cpp
>
> mkdir -p logs && env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
> nohup taskset -c 112 ./build/ecg1d_verify \
>   --N 4 --K 7 --run-step4 \
>   --T 6.283185307 --dt 1e-3 --dt-trace 1e-3 \
>   --rt-verbose --moment-form both \
>   > logs/realtime_N4_K7.log 2>&1 &
> ```
>
> Monitor: `tail -f logs/realtime_N4_K7.log`
> Check alive: `ps -o pid,etime,pcpu,cmd -p $(pgrep -f "ecg1d_verify --N 4")`
> Estimated runtime: ~2–4 hours (K²·N! ≈ 1200× the N=1 K=3 baseline of 5 s).
> When it finishes, ask me to make the figures.

## Worked example — sandbox-cwd variant

User: `/longtask realtime N=1 K=5 --no-wiener` while a default-Wiener N=1 K=5 job is already running on core 118.

Both jobs would write to the same `out/verify/step4_ecg_*_N1_K5.csv`. Sandbox the new one:

> Wrote `txt_task/realtime_N1_K5_no_wiener_cmds.txt`. Run:
>
> ```bash
> cd /home/gyqyan/zexuan/ecg1d_cpp
>
> mkdir -p runs/no_wiener_N1K5/out/verify logs
> ln -sf /home/gyqyan/zexuan/ecg1d_cpp/out/verify/step1_ecg_basis_N1_K5.csv \
>        runs/no_wiener_N1K5/out/verify/step1_ecg_basis_N1_K5.csv
>
> env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
> nohup taskset -c 113 bash -c '
>   cd /home/gyqyan/zexuan/ecg1d_cpp/runs/no_wiener_N1K5
>   /home/gyqyan/zexuan/ecg1d_cpp/build/ecg1d_verify \
>     --N 1 --K 5 --run-step4 \
>     --T 6.283185307 --dt 1e-3 --dt-trace 1e-3 \
>     --n-snap 11 --moment-form both --rt-verbose \
>     --no-wiener
> ' > logs/realtime_N1_K5_no_wiener.log 2>&1 &
> echo "PID=$!"
> ```
>
> Outputs land at `runs/no_wiener_N1K5/out/verify/step4_ecg_*_N1_K5.csv`, fully disjoint from the default-Wiener job's writes to `out/verify/`. Both jobs run concurrently on different cores and never collide.
