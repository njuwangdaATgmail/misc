set nocompatible
set backspace=indent,eol,start
if has('mouse')
  set mouse=a
endif
if has('gui_running')
  set guioptions=aiA
  set guifont=Consolas:h12
  set termguicolors
endif

set autoindent            "indent scheme: cindent, smartindent,...
set shiftwidth=2          "indent size
set tabstop=2             "tab size
set expandtab             "use space to expand tab

let fortran_do_enddo=1    "enable indent of do-enddo in fortran90

set showmatch     "match of bracket
set cursorline    "show a line at the cursor line

" plugin manager
call plug#begin('~/.vim/plugged')

  Plug 'lervag/vimtex'
  Plug 'altercation/vim-colors-solarized'
  "Plug 'arakashic/nvim-colors-solarized'
  Plug 'morhetz/gruvbox'
  Plug 'nightsense/snow'
  Plug 'andreypopp/vim-colors-plain'
  Plug 'Yggdroot/indentLine'

call plug#end()

" external pdf viewer used by vimtex
            let g:vimtex_view_general_viewer = 'SumatraPDF.exe'
            let g:vimtex_view_general_options
                \ = '-reuse-instance -forward-search @tex @line @pdf'
                \ . ' -inverse-search "gvim --servername ' . v:servername
                \ . ' --remote-send \"^<C-\^>^<C-n^>'
                \ . ':drop \%f^<CR^>:\%l^<CR^>:normal\! zzzv^<CR^>'
                \ . ':execute ''drop '' . fnameescape(''\%f'')^<CR^>'
                \ . ':\%l^<CR^>:normal\! zzzv^<CR^>'
                \ . ':call remote_foreground('''.v:servername.''')^<CR^>^<CR^>\""'

syntax enable
"colorscheme desert "default,desert,pablo,ron,solarized,torte,...
set background=dark
colorscheme solarized

" folding
set foldenable
set foldmethod=indent
let fortran_fold=1
set foldlevelstart=99

set textwidth=120
set colorcolumn=121
set statusline=%l,\ %v
set laststatus=2
