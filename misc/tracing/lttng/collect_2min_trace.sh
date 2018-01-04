#!/bin/bash

lttng create
lttng enable-event --userspace ubitrack:*
lttng start && sleep 120 && lttng stop
sudo chgrp -R $USER /home/$USER/lttng-traces/
sudo chmod -R g+w /home/$USER/lttng-traces/
