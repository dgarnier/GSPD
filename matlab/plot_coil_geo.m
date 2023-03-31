
function plot_coil_geo(tok_data_struct)

fcdata = tok_data_struct.fcdata;
fcnames = tok_data_struct.fcnames;

k = 17:22;
fcdata = fcdata(:,k);
fcnames = fcnames(k,:);

plot_polybox(fcdata(2,:),fcdata(1,:),fcdata(4,:),fcdata(3,:), [.5 .5 1],fcdata(5,:),fcdata(6,:))

for i = 1:length(fcnames)
text(fcdata(2,i), fcdata(1,i), fcnames(i,:), 'fontsize', 10, 'horizontalalignment', 'center')
end

end
