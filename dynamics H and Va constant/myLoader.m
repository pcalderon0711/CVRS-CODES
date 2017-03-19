function vec = myLoader(txt_file, type)

if type == 'i'
    dim = 12;
elseif type == 'p'
    dim = 36;
end

vec = zeros(dim,1);

ic_file = fopen(txt_file);
for i=1:dim
    vec(i) = str2double(fgetl(ic_file));
end