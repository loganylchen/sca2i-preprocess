#!/usr/bin/env bash
# Kill the given PID (and its process group) if free space on / drops below
# MIN_FREE_GB. Safety net for the per-cell full-genome reditools run, which
# accumulates large held-open temp during a long single process.
set -uo pipefail
TARGET_PID="${1:?usage: disk_watchdog.sh <pid> [min_free_gb]}"
MIN_FREE_GB="${2:-80}"
while kill -0 "$TARGET_PID" 2>/dev/null; do
    avail_gb=$(df -BG --output=avail / | tail -1 | tr -dc '0-9')
    if [ "${avail_gb:-999}" -lt "$MIN_FREE_GB" ]; then
        echo "$(date) WATCHDOG: only ${avail_gb}G free (< ${MIN_FREE_GB}G), killing $TARGET_PID" >&2
        kill -TERM -"$(ps -o pgid= "$TARGET_PID" | tr -d ' ')" 2>/dev/null || kill -TERM "$TARGET_PID"
        sleep 5
        kill -KILL "$TARGET_PID" 2>/dev/null
        exit 2
    fi
    sleep 60
done
echo "$(date) WATCHDOG: target $TARGET_PID exited; final free $(df -BG --output=avail / | tail -1 | tr -d ' ')" >&2
