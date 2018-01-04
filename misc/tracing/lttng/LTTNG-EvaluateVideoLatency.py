
# coding: utf-8

# In[81]:

import babeltrace
import pandas as pd
import numpy as np


# In[82]:

def parse_trace(session_name):
    tc1 = babeltrace.TraceCollection()
    tc1.add_trace("/home/netlabs/lttng-traces/{}/ust/uid/1000/64-bit".format(session_name), "ctf")
    msgs_start = {}
    msgs_latency = {}
    gpu_alloc = 0
    cpu_alloc = 0
    gpu_upload = 0
    gpu_download = 0
    for event in tc1.events:
        if event.get("port_field") == "VideoCapture":
            msgs_start[event.get("timestamp_field")] = event.timestamp
        elif event.get("port_field") == "TextureUpdated":
            if event.get("timestamp_field") in msgs_start.keys():
                msgs_latency[event.get("timestamp_field")] = event.timestamp - msgs_start[event.get("timestamp_field")]
        elif event.name == "ubitrack:vision_allocate_gpu":
            gpu_alloc += 1
        elif event.name == "ubitrack:vision_allocate_cpu":
            cpu_alloc += 1
        elif event.name == "ubitrack:vision_gpu_upload":
            gpu_upload += 1
        elif event.name == "ubitrack:vision_gpu_download":
            gpu_download += 1
    allocations = pd.DataFrame([{"cpu_alloc": cpu_alloc, "gpu_alloc": gpu_alloc,
                                "gpu_upload": gpu_upload, "gpu_download": gpu_download}])   
    latency = pd.DataFrame([msgs_latency], index=['duration_ms']).T/10**6 # to milliseconds
    return latency, allocations


# # FreenectGPU-Undistort-4xRotate90-BgImage

# In[83]:

latency_gpu, allocations_gpu = parse_trace("auto-20161209-193709")
print(latency_gpu.mean())
print(latency_gpu.describe())
print(allocations_gpu)


# # Freenect-Undistort-4xRotate90-BgImage

# In[87]:

# latency_cpu, allocations_cpu = parse_trace("auto-20161209-194045")
# print(latency_cpu.mean())
# print(latency_cpu.describe())
# print(allocations_cpu)

