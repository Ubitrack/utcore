#!/usr/sbin/dtrace -s

#pragma D option quiet

dtrace:::BEGIN
{
    printf("#Start Capturing Eventqueue Statistics, CTRL+C to show results\n");

}

ubitrack*:::eventqueue-dispatch-begin
{
    self->start = timestamp;
}

ubitrack*:::eventqueue-dispatch-end
{
    @times_component[arg0, copyinstr(arg2), copyinstr(arg3)] = sum(timestamp - self->start);
    @counts_component[arg0, copyinstr(arg2), copyinstr(arg3)] = count();
    @times_message[arg0, arg1] = sum(timestamp - self->start);
    @counts_message[arg0, arg1] = count();
}


END
{
        printf("#Time spent for pushing messages into component ports (microseconds)\n");
        printf(">times_component<");
        normalize(@times_component,1000);
        printa(@times_component);

        printf("#Number of messages processed by component ports\n");
        printf(">counts_component<");
        printa(@counts_component);

        printf("#Time spent this message priority by component ports (microseconds)\n");
        printf(">times_message<");
        normalize(@times_message,1000);
        printa(@times_message);

        printf("#Number of priorities processed by component ports\n");
        printf(">counts_message<");
        printa(@counts_message);

}
