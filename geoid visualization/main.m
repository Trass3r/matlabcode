gsigeome_ver5gsigeo20112000%  20.00000 120.00000 0.016667 0.025000 1801 1201 1%glamn = m_dataOffsetY - (m_nrows-1)/m_dataScaleY;%glomn = m_dataOffsetX;%dgla = 1.0 / m_dataScaleY;dgla = 0.016667;dglo = 0.025;%dglo = 1.0 / m_dataScaleX;nrows = 1801; % latitudencols = 1201; % longitudelatits = 20 : dgla : 20+dgla*(nrows-1);longits = 120 : dglo : 120+dglo*(ncols-1);% row-major -> column-major nastinessA = reshape(A, ncols, nrows)';Aold = reshape(Aold, ncols, nrows)';%A = flipud(A);%imagesc(longits, latits, A);imagesc([120 150], [20 50], A);pause(0.001); % jframe nastinessjframe = get(handle(gcf), 'JavaFrame');jframe.setMaximized(true); % maximize figurexlabel 'lon'; ylabel 'lat'; colorbar; axis xy; axis square;figure;imagesc([120 150], [20 50], A-Aold);pause(0.001); % jframe nastinessjframe = get(handle(gcf), 'JavaFrame');jframe.setMaximized(true); % maximize figurexlabel 'lon'; ylabel 'lat'; colorbar; axis xy; axis square;%caxis([0 100]);