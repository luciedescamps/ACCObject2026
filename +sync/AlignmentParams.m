classdef AlignmentParams
   
    properties
        % CHUNKSIZE number of consecutive pulses grouped as a chunk
        %
        % The alignment process works by dividing the first pulse train
        % into short "chunks" of equal length. Every chunk will contain a 
        % unique pattern of pulse intervals, so for every chunk we can find
        % the closest-matching chunk in the second pulse train.
        % 
        % A larger chunk size is more reliable, since it will reduce the 
        % probability of finding spurious matches. However, larger chunks
        % increase the requirement for pulses to be reliably detected in
        % both recordings, because any missing pulse in one pulse train
        % will invalidate the chunk enclosing the missing event. As an
        % example, consider a pulse train which failed to record 10% of the 
        % original pulses. Assuming the missed pulses are randomly
        % distributed, the probability that a chunk is valid is 0.9^10 =
        % 0.349. Increasing the chunk size to 30 yields a probability of
        % 0.9^30 = 0.0424.
        %
        % N.B. larger chunks take more time to align. Using more than
        chunkSize               (1,1)   double  {mustBeInteger, mustBeNonnegative} = 10
        
        % ALIGNMENTMETHOD specify the algorithm for matching chunks
        %
        % Two alternative algorithms are offered: "rg" (default) and 
        % "xcorr". "rg" matches a template chunk of pulse times to a chunk 
        % in the target pulse train, by finding the alignment which minimizes the squared
        % error between the pulse times of the template and target chunk. 
        % "xorr" uses chunks of inter-pulse intervals (IPIs) instead of 
        % pulse times, and finds the alignment which minimizes the mean
        % squared error in the IPIs.
        %
        % The "rg" method is generally more reliable, since it is sensitive
        % to pulse timing, which is cumulatively related to IPI. Changing 
        % an IPI in the middle of a chunk will shift the times of all 
        % following pulses in the chunk. The "rg" method would register
        % this cumulative effect on pulse times; conversely the "xcorr" 
        % method would only register the single IPI change, regardless of
        % its position in the chunk.
        alignmentMethod         (1,1)   string  {mustBeMember(alignmentMethod,["xcorr", "rg"])} = "rg"
        
        % ANCHOR pulse reference point used for alignment
        % 
        % Alignment can be defined in relation to the pulse rise times
        % "tRise", or fall times "tFall", or to the pulse mid-point, 
        % "tCenter" (default). "tCenter" may help to reduce alignment
        % errors if sampling quantization is a significant factor.
        anchor                  (1,1)   string  {mustBeMember(anchor, ["tRise", "tFall", "tCenter"])} = sync.BaseAligner.DEFAULT_ANCHOR
        
        % MINCHUNKSEPARATION chunk alignment quality criterion
        %
        % This parameter is a criterion for how well aligned a chunk must 
        % be for it to be kept. This quality check is necessary when the 
        % two pulse-train recordings do not represent an identical pulse
        % sequence. For example, there may be non-overlapping periods in
        % the two recordings, or there may be a small proportion of pulse
        % events which failed to registed in one recording. In these cases,
        % some chunks cannot be correctly matched between the two
        % recordings, and the best possible match will instead be spurious.
        % For this reason, we need to identify such mismatches and discard
        % them to prevent them from contaminating subsequent
        % synchronization operations.
        % 
        % Alignment quality for each chunk is measured by calculating the 
        % the mean squared error (MSE) for the alignment between the 
        % template chunk and the matched chunk. If the best alignment is 
        % correct, it should have a much smaller MSE than the second-best 
        % alignment. Conversely, if the correct alignment was not found and 
        % the best match is spurious, it should have a large MSE value 
        % which is similar to the second-best MSE value.
        %
        % If MSE1 is best MSE, and MSE2 is the second-best MSE, the 
        % separation S is defined as below: 
        % 
        %       S = log(MSE2) - log(MSE1) = log(MSE2/MSE1)
        % 
        % Therefore if S = 3.0, this means that MSE2 is larger than MSE1 by
        % a factor of exp(3.0) = 20.09.
        %
        % Alignments for which S is less than MINCHUNKSEPARATION will be
        % discarded.
        
        minChunkSeparation      (1,1)   double  {mustBeNonnegative} = 1
        
        % MINTOTALSEPARATION criterion for overall quality of alignment
        %
        % This parameter acts a quality check for the whole set of chunks.
        % The single-chunk separation values are averaged to calculate the
        % total separation. The value of MINTOTALSEPARATION is the lowest
        % permissible value for the total separation. A total separation of
        % less than this will result in an error being thrown by the
        % Aligner.align() method.
        minTotalSeparation      (1,1)   double  {mustBeNonnegative} = 2.5


        
        verbose                 (1,1)   logical = false

        truncate                (1,1)   logical = true
    end
    
    methods
        function ok = checkConsistency(self, params)
            % Check that values are consistent with another parameter set
            fds = fieldnames(self);
            ok = true;
            for f = 1:numel(fds)
                fd = fds{f};
                v1 = self.(fd);
                v2 = params.(fd);
                if ~isequal(v1, v2)
                    warning('AlignmentParams:checkConsistency', ...
                        'Inconsistent value for parameter "%s" (%s vs %s)', ...
                        fd, string(v1), string(v2));
                    ok = false;
                end
            end
        end
    end
    
end