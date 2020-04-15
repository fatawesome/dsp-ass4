clf()
clear()

function sv = shift(v,n)
if (size(v,'r')<>1 & size(v,'c')<>1) then
    error("1st argument must be column vector or row vector")
end
if size(v,'r')<>1 then
    v = shift(v',n)
    sv = v'
else
    n = modulo(n,size(v,'c'))
    sv = v($-n+1:$)
    sv = [sv v(1:$-n)]
end
endfunction

[irc, Fs_irc, _] = wavread("./fir.wav");
[poetry, Fs_poetry, _] = wavread("./poetry.wav");
irc = irc(1, :)
irc = irc(1, 182000:length(irc));
poetry = poetry(1, :)
//plot(abs(fft(poetry)))

//irc = shift(irc', (length(irc)-1)/2)'

[db, phi] = dbphi(fft(irc));

inverse_filter = ifft(conj(fft(irc))./abs(fft(irc)))
//plot(inverse_filter)
inverse_filter = shift(inverse_filter', (length(inverse_filter)-1)/2)'
inverse_filter = inverse_filter .* window('kr', length(inverse_filter), 8)

[db, phi] = dbphi(fft(inverse_filter));
//plot(phi)


modified = convol(poetry, irc)
cleaned = convol(modified, inverse_filter)
plot(convol(irc, inverse_filter))

savewave('my-irc_modified.wav', modified, 44100)
savewave('my-irc_cleaned.wav', cleaned, 44100)
