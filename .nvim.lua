vim.keymap.set('n', '<localleader>t', function()
	-- run("clear && make test_project && clear && ./test_project")
	require("snacks").terminal.get("make test_project && clear && ./test_project", { interactive = false }):show()
end)

vim.keymap.set('n', '<localleader>m', function()
	-- run("clear && make app && clear && ./app")
	require("snacks").terminal.get("make app && clear && ./app")
end)
