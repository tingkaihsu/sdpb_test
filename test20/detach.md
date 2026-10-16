```bash
nohup bash ./run_sdpb.sh > sdpb.log 2>&1 < /dev/null &

echo $!
```
Save the background job PID.

Detach the shell from the job
```bash
disown
```

Reconnection
```bash
tail -f sdpb.log
```
or
```bash
ps -p <PID>
```

PID: `1741746`
