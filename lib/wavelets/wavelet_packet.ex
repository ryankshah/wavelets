defmodule Wavelets.WaveletPacket do
  @moduledoc """
  Implementation of Wavelet Packet decomposition and reconstruction.
  """

  alias Wavelets.DWT
  alias Wavelets.Filter
  alias Wavelets.IDWT

  @doc """
  Performs 1D wavelet packet decomposition
  """
  @spec decompose_1d(
          list(number),
          Filter.t(),
          non_neg_integer(),
          Wavelets.precision()
        ) ::
          map()
  def decompose_1d(signal, filter, levels, _precision \\ :double) do
    initial_tree = %{{0, 0} => signal}
    do_decompose(initial_tree, filter, 0, levels)
  end

  defp do_decompose(tree, _filter, current_level, max_levels)
       when current_level >= max_levels do
    tree
  end

  defp do_decompose(tree, filter, current_level, max_levels) do
    new_nodes =
      for j <- 0..(:math.pow(2, current_level) |> trunc |> Kernel.-(1)),
          signal when is_list(signal) <- [Map.get(tree, {current_level, j})],
          {approx, details} = DWT.forward_1d(signal, filter) do
        %{
          {current_level + 1, 2 * j} => approx,
          {current_level + 1, 2 * j + 1} => details
        }
      end
      |> Enum.reduce(%{}, &Map.merge(&2, &1))

    tree = Map.merge(tree, new_nodes)
    do_decompose(tree, filter, current_level + 1, max_levels)
  end

  @doc """
  Reconstructs signal from wavelet packet decomposition
  """
  @spec reconstruct_1d(map(), Filter.t(), Wavelets.precision()) :: list(number)
  def reconstruct_1d(tree, filter, _precision \\ :double) do
    max_level =
      tree
      |> Map.keys()
      |> Enum.map(fn {level, _} -> level end)
      |> Enum.max()

    reconstruct_level(tree, filter, max_level)
  end

  defp reconstruct_level(tree, _filter, 0), do: Map.get(tree, {0, 0})

  defp reconstruct_level(tree, filter, level) do
    new_nodes =
      0..(:math.pow(2, level - 1) |> trunc |> Kernel.-(1))
      |> Enum.map(fn j ->
        approx = Map.get(tree, {level, 2 * j})
        details = Map.get(tree, {level, 2 * j + 1})

        if approx && details do
          node = IDWT.inverse_1d(approx, details, filter)
          {{level - 1, j}, node}
        end
      end)
      |> Enum.reject(&is_nil/1)
      |> Map.new()

    tree = Map.merge(tree, new_nodes)
    reconstruct_level(tree, filter, level - 1)
  end
end
