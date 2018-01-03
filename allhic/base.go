/**
 * Filename: /Users/bao/code/allhic/allhic/base.go
 * Path: /Users/bao/code/allhic/allhic
 * Created Date: Tuesday, January 2nd 2018, 8:07:22 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import (
	"os"
	"path"
	"strings"

	logging "github.com/op/go-logging"
)

var log = logging.MustGetLogger("allhic")
var format = logging.MustStringFormatter(
	`%{color}%{time:15:04:05} %{shortfunc} ▶ %{level:.4s} %{color:reset} %{message}`,
)

// Backend is the default stderr output
var Backend = logging.NewLogBackend(os.Stderr, "", 0)

// BackendFormatter contains the fancy debug formatter
var BackendFormatter = logging.NewBackendFormatter(Backend, format)

// RemoveExt returns the substring minus the extension
func RemoveExt(filename string) string {
	return strings.TrimSuffix(filename, path.Ext(filename))
}
