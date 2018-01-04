from etw import TraceEventSource, EventConsumer, EventHandler
from etw.descriptors import image
import exceptions
import os
import sys
import numpy as np

from etw.descriptors import uteventqueue

class NodeLatencyConsumer(EventConsumer):

    def __init__(self):
        self.eventqueue_stats = {}

    # @EventHandler(uteventqueue.Event.DispatchBegin)
    # def OnDispatchBegin(self, event_data):
    #   print event_data.time_stamp, event_data.EventDomain, event_data.Priority, event_data.ComponentName, event_data.PortName

    @EventHandler(uteventqueue.Event.DispatchEnd)
    def OnDispatchEnd(self, event_data):
        key = (event_data.EventDomain, event_data.ComponentName, event_data.PortName)
        samples = self.eventqueue_stats.setdefault(key, [])
        samples.append(np.asarray((event_data.time_stamp, event_data.Priority, event_data.Duration)))



class AllocationConsumer(EventConsumer):

    def __init__(self):
        self.eventqueue_stats = {}

    @EventHandler(uteventqueue.Event.VisionAllocateGpu)
    def OnVisionAllocateGpu(self, event_data):
        key = "allocate_gpu"
        samples = self.eventqueue_stats.setdefault(key, 0)
        self.eventqueue_stats[key] = samples + 1

    @EventHandler(uteventqueue.Event.VisionAllocateCpu)
    def OnVisionAllocateCpu(self, event_data):
        key = "allocate_cpu"
        samples = self.eventqueue_stats.setdefault(key, 0)
        self.eventqueue_stats[key] = samples + 1


    @EventHandler(uteventqueue.Event.VisionGpuUpload)
    def OnVisionGpuUpload(self, event_data):
        key = "upload_gpu"
        samples = self.eventqueue_stats.setdefault(key, 0)
        self.eventqueue_stats[key] = samples + 1

    @EventHandler(uteventqueue.Event.VisionGpuDownload)
    def OnVisionGpuDownload(self, event_data):
        key = "download_gpu"
        samples = self.eventqueue_stats.setdefault(key, 0)
        self.eventqueue_stats[key] = samples + 1


class MessageLatencyConsumer(EventConsumer):

    def __init__(self):
        self.eventqueue_stats = {}
        self.create_ts = {}

    @EventHandler(uteventqueue.Event.CreateMeasurement)
    def OnCreateMeasurement(self, event_data):
        key = event_data.Priority
        # will not work with complex graphs!!!
        self.create_ts[key] = event_data.raw_time_stamp
        
    @EventHandler(uteventqueue.Event.ReceiveMeasurement)
    def OnReceiveMeasurement(self, event_data):
        key = event_data.Priority
        cts = self.create_ts.get(key, event_data.raw_time_stamp)
        dt = event_data.raw_time_stamp - cts
        self.eventqueue_stats[key] = dt


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "usage: %s <trace>.etl" % sys.argv[0]
        exit(0)

    filename = sys.argv[1]
    ts = NodeLatencyConsumer()
    ac = AllocationConsumer()
    mc = MessageLatencyConsumer()
    tes = TraceEventSource([ts, ac, mc], True)

    print "open tracefile"
    tes.OpenFileSession(filename)

    print"parsing entries ..."
    tes.Consume()

    print "average time (ms)"
    for k in sorted(ts.eventqueue_stats.keys()):
        values = ts.eventqueue_stats[k]
        data = np.vstack(values)
        ed, cn, pn = k
        print "%s:%s-%s (%d) = %0.5f" % (ed, cn, pn, len(values), data[:,2].mean())

    print "allocations"
    for k in sorted(ac.eventqueue_stats.keys()):
        value = ac.eventqueue_stats[k]
        print "%s = %d" % (k, value)

    data = np.array(mc.eventqueue_stats.values())
    print "message_latency (ms): %0.5f" % (data.mean() * (1.0/10000.0))

