function vec = myLoader(txt_file, type)

if type == 'i'
    dim = 14;
elseif type == 'p'
    dim = 34;
end

vec = zeros(dim,1);

ic_file = fopen(txt_file);
for i=1:dim
    vec(i) = str2double(fgetl(ic_file));
end