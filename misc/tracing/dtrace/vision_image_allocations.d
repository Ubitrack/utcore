#!/usr/sbin/dtrace -s

#pragma D option quiet

dtrace:::BEGIN
{
    printf("#Start Capturing Image Allocation Statistics, CTRL+C to show results\n");
}

ubitrack*:::vision-allocate-cpu
{
    @counts_alloc_cpu[arg0] = count();
}

ubitrack*:::vision-allocate-gpu
{
    @counts_alloc_gpu[arg0] = count();
}

END
{
        printf("#Count allocations\n");
        printf(">CPU<");
        printa(@counts_alloc_cpu);
        printf(">GPU<");
        printa(@counts_alloc_gpu);
}
