#!/usr/bin/perl -w

use threads;
use Thread::Queue;

my @stufftodo = (1..1000);

my $THREADS= 1; # Number of threads

my $workq = Thread::Queue->new(); # Work to do

$workq->enqueue(@stufftodo); # Queue up some work to do
$workq->enqueue("EXIT") for(1..$THREADS); # And tell them when
                                          # they're done

threads->create("Handle_Work") for(1..$THREADS); # Spawn our workers

# Process returns while the the threads are running.
# Alternatively, if we just want to wait until they're all done:
# sleep 10 while threads->list(threads::running);
while(threads->list(threads::running)){	
   # Do something with the data
   # Salt to taste
}
# When we get here, there are no more running threads.
# At this point we may want to take one more run through the
# return queue, or do whatever makes sense.

sub Handle_Work {
  while(my $todo=$workq->dequeue()) {
    last if $todo eq 'EXIT'; # All done
	print $todo."\n";
  }

  # We're all done with this thread
  threads->detach;
}
