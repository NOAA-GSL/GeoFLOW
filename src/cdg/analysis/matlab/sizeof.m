function bytes = sizeof(x)
    w = whos('x');
    bytes = w.bytes / numel(x);
end
