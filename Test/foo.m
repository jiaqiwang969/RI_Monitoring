function out = foo(in)
    validateattributes(in,{'numeric'},{'nonempty'}); %Not now
    % Returns zero
    out = zeros(size(in),'like',in);
end